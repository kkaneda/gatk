package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

/** Class for storing a segment that was generated by a model (e.g. RCS, ACS), though we do not necessarily know ahead of time
 *
 * @author lichtens &lt;lichtens@broadinstitute.org&gt;
 */
public class ModeledSegment implements Locatable {
    protected SimpleInterval simpleInterval;

    /** Segment Mean in log copy ratio space as determined by whatever model generated this segment.
     */
    private double segmentMean;

    private String call;

    public ModeledSegment(final SimpleInterval interval, final String call, final double segmentMean) {
        this.simpleInterval = Utils.nonNull(interval, "The input interval cannot be null");
        this.call = Utils.nonNull(call, "The input call cannot be null.  Use empty string, instead \"\" ");
        this.segmentMean = segmentMean;
    }

    public double getSegmentMean() {
        return segmentMean;
    }

    public void setSegmentMean(final double segmentMean) {
        this.segmentMean = segmentMean;
    }

    public double getSegmentMeanInCRSpace() {
        return Math.pow(2, segmentMean);
    }

    public void setSegmentMeanInCRSpace(final double segmentMeanInCRSpace) {
        this.segmentMean = Math.log(segmentMeanInCRSpace)/Math.log(2);
    }

    @Override
    public String getContig() {return simpleInterval.getContig(); }

    @Override
    public int getStart() {return simpleInterval.getStart(); }

    @Override
    public int getEnd() {return simpleInterval.getEnd(); }

    public SimpleInterval getSimpleInterval() {
        return simpleInterval;
    }

    public void setSimpleInterval(final SimpleInterval simpleInterval) {
        this.simpleInterval = Utils.nonNull(simpleInterval, "The input interval cannot be null");
    }

    /**
     * Returns the call.  Returns null for uncalled segments.
     *
     * @return never {@code null}
     */
    public String getCall() {
        return this.call;
    }

    /**
     * Sets the call.
     */
    public void setCall(final String call) {
        this.call = Utils.nonNull(call, "The input call cannot be null.  Use empty string, instead \"\" ");
    }
}