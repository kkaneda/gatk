package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

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

    @Override
    public boolean requiresReads()
    {
        return true;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {

        final Broadcast<Map<String, Integer>> broadcastContigNameToIDMap =
                ctx.broadcast(ReadMetadata.buildContigNameToIDMap(getHeaderForReads()));

        JavaRDD<GATKRead> reads =
                getUnfilteredReads()
                    .filter(FindSmallIndelRegions::isTestableRead);

        final List<IntHistogram> histograms =
                reads
                    .mapPartitions(readItr -> {
                        final IntHistogram histogram = new IntHistogram(MAX_TRACKED_VALUE);
                        while ( readItr.hasNext() ) {
                            histogram.addObservation(Math.abs(readItr.next().getFragmentLength()));
                        }
                        return Collections.singletonList(histogram).iterator();
                    })
                    .collect();

        final IntHistogram sumHistogram = new IntHistogram(MAX_TRACKED_VALUE);
        for ( final IntHistogram histogram : histograms ) {
            sumHistogram.addObservations(histogram);
        }
        final Broadcast<IntHistogram.CDF> broadcastCDF = ctx.broadcast(sumHistogram.trim().getCDF());

        final List<BreakpointEvidence> evidenceList =
                reads
                    .mapPartitions(readItr -> {
                        if ( !readItr.hasNext() ) return Collections.emptyIterator();
                        final List<BreakpointEvidence> evList = new ArrayList<>();
                        final IntHistogram.CDF cdf = broadcastCDF.value();
                        final Map<String, Integer> contigNameToIDMap = broadcastContigNameToIDMap.value();
                        final IntHistogram[] histoPair = new IntHistogram[2];
                        histoPair[0] = cdf.createEmptyHistogram();
                        histoPair[1] = cdf.createEmptyHistogram();
                        GATKRead read = readItr.next();
                        String curContig = read.getContig();
                        int curEnd = read.getStart() + BLOCK_SIZE;
                        int fillIdx = 0;
                        do {
                            if ( !read.getContig().equals(curContig) || read.getStart() >= curEnd ) {
                                final int oldIdx = fillIdx ^ 1; // flip the lowest bit
                                histoPair[oldIdx].addObservations(histoPair[fillIdx]);
                                if ( cdf.isDifferent(histoPair[oldIdx]) ) {
                                    evList.add(createEvidence(contigNameToIDMap.get(curContig), curEnd));
                                }
                                histoPair[oldIdx].clear();
                                fillIdx = oldIdx;
                                if ( read.getContig().equals(curContig) ) {
                                    curEnd += BLOCK_SIZE;
                                } else {
                                    curContig = read.getContig();
                                    curEnd = read.getStart() + BLOCK_SIZE;
                                }
                            }
                            histoPair[fillIdx].addObservation(Math.abs(read.getFragmentLength()));
                        } while ( (read = readItr.hasNext() ? readItr.next() : null) != null );

                        final int oldIdx = fillIdx ^ 1; // flip the lowest bit
                        histoPair[oldIdx].addObservations(histoPair[fillIdx]);
                        if ( cdf.isDifferent(histoPair[oldIdx]) ) {
                            evList.add(createEvidence(contigNameToIDMap.get(curContig), curEnd));
                        }
                        return evList.iterator();
                    })
                    .collect();
        for ( final BreakpointEvidence evidence : evidenceList ) {
            System.out.println(evidence);
        }
    }

    private static boolean isTestableRead( final GATKRead read ) {
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
