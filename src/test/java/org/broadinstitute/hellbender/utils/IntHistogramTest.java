package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.Random;

public class IntHistogramTest extends BaseTest {
    private static final int MAX_TRACKED_VALUE = 1000;
    private static final float SIGNIFICANCE = .05f;
    private static Random random;

    @BeforeMethod
    public void reseedRandom() {
        random = new Random(47L);
    }

    @Test
    void testTrim() {
        final IntHistogram intHistogram = new IntHistogram(8);
        intHistogram.addObservation(0);
        intHistogram.addObservations(1, 3L);
        intHistogram.addObservations(2, 5L);
        intHistogram.addObservations(3, 3L);
        intHistogram.addObservation(4);
        intHistogram.addObservation(10);
        Assert.assertEquals(intHistogram.trim().getMaximumTrackedValue(), 3);
    }

    @Test
    void testTrimTo0() {
        final IntHistogram intHistogram = new IntHistogram(2);
        intHistogram.addObservations(0, 10L);
        intHistogram.addObservation(3);
        Assert.assertEquals(intHistogram.trim().getMaximumTrackedValue(), 0);
    }

    @Test
    void testTrimEarlyBailout() {
        final IntHistogram intHistogram = new IntHistogram(0);
        intHistogram.addObservations(0, 10L);
        intHistogram.addObservations(1, 10L);
        Assert.assertTrue(intHistogram.trim() == intHistogram);
    }

    // TODO: figure out what the data should really look like
    @Test
    public void testSmallInsertion() {
        final IntHistogram.CDF cdf = genNormalSample(480, 25, 10000000).getCDF();
        final IntHistogram sample = genNormalSample(480, 25, 25);
        for ( int idx = 75; idx != MAX_TRACKED_VALUE; ++idx )
            sample.addObservations(idx-50, sample.getNObservations(idx));
        Assert.assertTrue(cdf.isDifferent(sample, .01f));
    }

    @Test
    public void testNormalDistribution() {
        final IntHistogram.CDF cdf = genNormalSample(480, 25, 10000).getCDF();
        Assert.assertFalse(cdf.isDifferent(genNormalSample(480, 25, 100), SIGNIFICANCE));
        // these numbers were jiggled by hand until the difference was barely detectable
        Assert.assertTrue(cdf.isDifferent(genNormalSample(489, 25, 100), SIGNIFICANCE));
        Assert.assertTrue(cdf.isDifferent(genNormalSample(473, 25, 100), SIGNIFICANCE));
        Assert.assertTrue(cdf.isDifferent(genNormalSample(480, 15, 100), SIGNIFICANCE));
        Assert.assertTrue(cdf.isDifferent(genNormalSample(480, 42, 100), SIGNIFICANCE));
    }

    @Test
    public void testLogNormalDistribution() {
        final IntHistogram.CDF cdf = genNormalSample(480, 25, 10000).getCDF();
        Assert.assertFalse(cdf.isDifferent(genLogNormalSample(480, 25, 100), SIGNIFICANCE));
        Assert.assertTrue(cdf.isDifferent(genLogNormalSample(490, 25, 100), SIGNIFICANCE));
        Assert.assertTrue(cdf.isDifferent(genLogNormalSample(473, 25, 100), SIGNIFICANCE));
        Assert.assertTrue(cdf.isDifferent(genLogNormalSample(480, 15, 100), SIGNIFICANCE));
        Assert.assertTrue(cdf.isDifferent(genLogNormalSample(480, 42, 100), SIGNIFICANCE));
    }

    @Test
    public void testBiModalDistribution() {
        final IntHistogram.CDF cdf = genNormalSample(480, 25, 10000).getCDF();
        final IntHistogram sample1 = genNormalSample(480, 25, 50);
        final IntHistogram sample2 = genNormalSample(500, 25, 50);
        sample1.addObservations(sample2);
        Assert.assertTrue(cdf.isDifferent(sample1, SIGNIFICANCE));
    }

    private static IntHistogram genNormalSample( final int mean, final int stdDev, final int nSamples ) {
        final IntHistogram histogram = new IntHistogram(MAX_TRACKED_VALUE);
        for ( int sample = 0; sample != nSamples; ++sample ) {
            histogram.addObservation((int)Math.round(mean + stdDev * random.nextGaussian()));
        }
        return histogram;
    }

    private static IntHistogram genLogNormalSample( final int mean, final int stdDev, final int nSamples ) {
        final double phi = Math.sqrt((double)stdDev * stdDev + (double)mean * mean);
        final double mu = Math.log((double)mean * mean / phi);
        final double sigma = Math.sqrt(Math.log(phi * phi / mean / mean));
        final IntHistogram histogram = new IntHistogram(MAX_TRACKED_VALUE);
        for ( int sample = 0; sample != nSamples; ++sample ) {
            histogram.addObservation((int)Math.round(Math.exp(mu + sigma * random.nextGaussian())));
        }
        return histogram;
    }
}