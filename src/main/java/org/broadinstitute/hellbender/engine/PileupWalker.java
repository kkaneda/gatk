package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.downsampling.PositionalDownsampler;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * An AssemblyRegionWalker is a tool that processes an entire region of reads at a time, each marked as either "active"
 * (containing possible variation) or "inactive" (not likely to contain actual variation). Tool authors must implement
 * {@link #apply} to process each region. Authors must also implement methods providing default values for the
 * various traversal parameters.
 *
 * {@link #apply} will be called once for each active AND inactive region, and it is up to the implementation how to
 * handle/process active vs. inactive regions.
 *
 * Internally, the reads are loaded in chunks called read shards, which are then subdivided into active/inactive regions
 * for processing by the tool implementation. Read shards should typically be much larger than the maximum assembly
 * region size to achieve good performance, and should have sufficient padding on either end to avoid boundary artifacts
 * for events near the shard boundaries.
 *
 * Read shards exist mainly as a proof-of-concept that we can shard the reads without introducing calling artifacts,
 * which will be important for the Spark equivalent of this traversal.
 */
public abstract class PileupWalker extends GATKTool {

    @Argument(fullName="readShardSize", shortName="readShardSize", doc = "Maximum size of each read shard, in bases. For good performance, this should be much larger than the maximum assembly region size.", optional = true)
    protected int readShardSize = defaultReadShardSize();

    @Argument(fullName="readShardPadding", shortName="readShardPadding", doc = "Each read shard has this many bases of extra context on each side. Read shards must have as much or more padding than assembly regions.", optional = true)
    protected int readShardPadding = defaultReadShardPadding();

    @Argument(fullName = "maxReadsPerAlignmentStart", shortName = "maxReadsPerAlignmentStart", doc = "Maximum number of reads to retain per alignment start position. Reads above this threshold will be downsampled. Set to 0 to disable.", optional = true)
    protected int maxReadsPerAlignmentStart = defaultMaxReadsPerAlignmentStart();

    /**
     * @return Default value for the {@link #readShardSize} parameter, if none is provided on the command line
     */
    protected abstract int defaultReadShardSize();

    /**
     * @return Default value for the {@link #readShardPadding} parameter, if none is provided on the command line
     */
    protected abstract int defaultReadShardPadding();

    /**
     * @return Default value for the {@link #maxReadsPerAlignmentStart} parameter, if none is provided on the command line
     */
    protected abstract int defaultMaxReadsPerAlignmentStart();

    @Override
    public final boolean requiresReads() { return true; }

    @Override
    public final boolean requiresReference() { return true; }

    @Override
    public String getProgressMeterRecordLabel() { return "regions"; }
    
    private List<LocalReadShard> readShards;

    /**
     * Initialize data sources for traversal.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalStart() instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();

        if ( readShardSize <= 0 ) {
            throw new CommandLineException.BadArgumentValue("read shard size must be > 0");
        }

        if ( readShardPadding < 0 ) {
            throw new CommandLineException.BadArgumentValue("read shard padding must be >= 0");
        }

        final List<SimpleInterval> intervals = hasIntervals() ? intervalsForTraversal : IntervalUtils.getAllIntervalsForReference(getHeaderForReads().getSequenceDictionary());
        readShards = makeReadShards(intervals);
    }

    /**
     * Shard our intervals for traversal into ReadShards using the {@link #readShardSize} and {@link #readShardPadding} arguments
     *
     * @param intervals unmodified intervals for traversal
     * @return List of {@link LocalReadShard} objects, sharded and padded as necessary
     */
    private List<LocalReadShard> makeReadShards(final List<SimpleInterval> intervals ) {
        final List<LocalReadShard> shards = new ArrayList<>();

        for ( final SimpleInterval interval : intervals ) {
            shards.addAll(LocalReadShard.divideIntervalIntoShards(interval, readShardSize, readShardPadding, reads, getHeaderForReads().getSequenceDictionary()));
        }

        return shards;
    }

    /**
     * Returns the default list of CommandLineReadFilters that are used for this tool. The filters
     * returned by this method are subject to selective enabling/disabling and customization by the
     * user via the command line. The default implementation uses the {@link WellformedReadFilter}
     * filter with all default options, as well as the {@link ReadFilterLibrary.MappedReadFilter}.
     * Subclasses can override to provide alternative filters.
     *
     * Note: this method is called before command line parsing begins, and thus before a SAMFileHeader is
     * available through {link #getHeaderForReads}.
     *
     * @return List of default filter instances to be applied for this tool.
     */
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> defaultFilters = new ArrayList<>(2);
        defaultFilters.add(new WellformedReadFilter());
        defaultFilters.add(new ReadFilterLibrary.MappedReadFilter());
        return defaultFilters;
    }

    @Override
    public final void traverse() {

        CountingReadFilter countedFilter = makeReadFilter();

        // Since we're processing regions rather than individual reads, tell the progress
        // meter to check the time more frequently (every 10 regions instead of every 1000 regions).
        progressMeter.setRecordsBetweenTimeChecks(10L);

        for ( final LocalReadShard readShard : readShards ) {
            // Since reads in each shard are lazily fetched, we need to pass the filter to the window
            // instead of filtering the reads directly here
            readShard.setReadFilter(countedFilter);
            readShard.setDownsampler(maxReadsPerAlignmentStart > 0 ? new PositionalDownsampler(maxReadsPerAlignmentStart, getHeaderForReads()) : null);

            final SAMFileHeader readsHeader = getHeaderForReads();
            final ReferenceContext referenceContext = new ReferenceContext(reference, readShard.getPaddedInterval());
            final FeatureContext featureContext = new FeatureContext(features, readShard.getPaddedInterval());
            final List<GATKRead> windowReads = Utils.stream(readShard).collect(Collectors.toList());
            final LocusIteratorByState locusIterator = new LocusIteratorByState(windowReads.iterator(), DownsamplingMethod.NONE,
                    false, ReadUtils.getSamplesFromHeader(readsHeader), readsHeader, false);
            locusIterator.forEachRemaining(alignmentContext -> apply(alignmentContext, referenceContext, featureContext));
        }

        logger.info(countedFilter.getSummaryLine());
    }

    /**
     * Shutdown data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalDone() instead.
     */
    @Override
    protected final void onShutdown() {
        // Overridden only to make final so that concrete tool implementations don't override
        super.onShutdown();
    }

    /**
     * Process an individual ReadPileup. Must be implemented by tool authors.
     *
     * @param alignmentContext ReadPileup plus location to process
     * @param referenceContext reference data overlapping the full extended span of the assembly region
     * @param featureContext features overlapping the full extended span of the assembly region
     */
    public abstract void apply( final AlignmentContext alignmentContext, final ReferenceContext referenceContext, final FeatureContext featureContext );
}
