package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;

@CommandLineProgramProperties(summary="Find regions likely to contain small indels due to fragment length anomalies.",
        oneLineSummary="Find regions containing small indels.",
        programGroup=StructuralVariationSparkProgramGroup.class)
public final class FindSmallIndelRegions extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private static final int BLOCK_SIZE = 250;
    private static final int MAX_TRACKED_VALUE = 2000;
    private static final int EVIDENCE_WEIGHT = 10;
    private static final int MIN_MAPQ = 20;
    private static final int MIN_MATCH_LEN = 45;
    private static final int ALLOWED_SHORT_FRAGMENT_OVERHANG = 10;
    private static final int EVIDENCE_SIZE_GUESS = 1000;

    @Override
    public boolean requiresReads()
    {
        return true;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {

        final JavaRDD<GATKRead> allReads = getUnfilteredReads();
        final JavaRDD<GATKRead> mappedReads =
                allReads.filter(read -> !read.isDuplicate() && !read.failsVendorQualityCheck() && !read.isUnmapped());

        final ReadMetadata readMetadata =
                new ReadMetadata(Collections.emptySet(), getHeaderForReads(), MAX_TRACKED_VALUE, mappedReads);
        final Broadcast<ReadMetadata> broadcastReadMetadata = ctx.broadcast(readMetadata);

        final JavaRDD<GATKRead> testableReads = allReads.filter(FindSmallIndelRegions::isTestableRead);

        final List<BreakpointEvidence> evidenceList =
            testableReads
                .mapPartitions(readItr -> gatherEvidence(readItr,broadcastReadMetadata.value()))
                .collect();

        for ( final BreakpointEvidence evidence : evidenceList ) {
            System.out.println(evidence);
        }
    }

    public static Iterator<BreakpointEvidence> gatherEvidence( final Iterator<GATKRead> readItr,
                                                               final ReadMetadata readMetadata ) {
        if ( !readItr.hasNext() ) return Collections.emptyIterator();

        final List<BreakpointEvidence> evList = new ArrayList<>(EVIDENCE_SIZE_GUESS);
        final Map<String, IntHistogram[]> libraryToHistoPairMap =
                new HashMap<>(SVUtils.hashMapCapacity(readMetadata.getAllLibraryStatistics().size()));
        GATKRead read = readItr.next();
        String curContig = read.getContig();
        int curEnd = read.getStart() + BLOCK_SIZE;
        int fillIdx = 0;
        do {
            if ( !read.getContig().equals(curContig) || read.getStart() >= curEnd ) {
                final int oldIdx = fillIdx ^ 1; // flip the lowest bit
                for ( final Map.Entry<String,IntHistogram[]> entry : libraryToHistoPairMap.entrySet() ) {
                    final IntHistogram[] histoPair = entry.getValue();
                    histoPair[oldIdx].addObservations(histoPair[fillIdx]);
                    if ( readMetadata.getLibraryStatistics(entry.getKey()).isDifferent(histoPair[oldIdx]) ) {
                        evList.add(createEvidence(readMetadata.getContigID(curContig), curEnd));
                    }
                    histoPair[oldIdx].clear();
                }
                fillIdx = oldIdx;
                if ( read.getContig().equals(curContig) ) {
                    curEnd += BLOCK_SIZE;
                } else {
                    curContig = read.getContig();
                    curEnd = read.getStart() + BLOCK_SIZE;
                }
            }
            final IntHistogram[] histoPair =
                libraryToHistoPairMap.computeIfAbsent(readMetadata.getLibraryName(read.getReadGroup()),
                            libName -> {
                               final IntHistogram[] pair = new IntHistogram[2];
                               final IntHistogram.CDF cdf = readMetadata.getLibraryStatistics(libName);
                                pair[0] = cdf.createEmptyHistogram();
                                pair[1] = cdf.createEmptyHistogram();
                               return pair; });
            histoPair[fillIdx].addObservation(Math.abs(read.getFragmentLength()));
        } while ( (read = readItr.hasNext() ? readItr.next() : null) != null );

        final int oldIdx = fillIdx ^ 1; // flip the lowest bit
        for ( final Map.Entry<String,IntHistogram[]> entry : libraryToHistoPairMap.entrySet() ) {
            final IntHistogram[] histoPair = entry.getValue();
            histoPair[oldIdx].addObservations(histoPair[fillIdx]);
            if ( readMetadata.getLibraryStatistics(entry.getKey()).isDifferent(histoPair[oldIdx]) ) {
                evList.add(createEvidence(readMetadata.getContigID(curContig), curEnd));
            }
        }
        return evList.iterator();
    }

    public static boolean isTestableRead( final GATKRead read ) {
        return read.isFirstOfPair() &&
                !read.isUnmapped() &&
                !read.mateIsUnmapped() &&
                read.isReverseStrand() != read.mateIsReverseStrand() &&
                !read.isDuplicate() &&
                !read.failsVendorQualityCheck() &&
                !read.isSecondaryAlignment() &&
                !read.isSupplementaryAlignment() &&
                read.getMappingQuality() >= MIN_MAPQ &&
                SVUtils.matchLen(read.getCigar()) >= MIN_MATCH_LEN &&
                read.getContig().equals(read.getMateContig()) &&
                (read.isReverseStrand() ?
                    read.getStart() + ALLOWED_SHORT_FRAGMENT_OVERHANG < read.getMateStart() :
                    read.getStart() - ALLOWED_SHORT_FRAGMENT_OVERHANG > read.getMateStart());
    }

    private static BreakpointEvidence createEvidence( final int contigId, final int end ) {
        int start = end - 2*BLOCK_SIZE;
        if ( start < 1 ) start = 1;
        final SVInterval interval = new SVInterval(contigId, start, end);
        return new BreakpointEvidence.TemplateSizeAnomaly(interval,EVIDENCE_WEIGHT);
    }
}
