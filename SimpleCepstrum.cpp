/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */
/*
    This file is Copyright (c) 2012 Chris Cannam
  
    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use, copy,
    modify, merge, publish, distribute, sublicense, and/or sell copies
    of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
    ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
    CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
    WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "SimpleCepstrum.h"

#include "vamp-sdk/FFT.h"

#include <vector>
#include <algorithm>

#include <cstdio>
#include <cmath>
#include <complex>

#if ( VAMP_SDK_MAJOR_VERSION < 2 || ( VAMP_SDK_MAJOR_VERSION == 2 && VAMP_SDK_MINOR_VERSION < 4 ) )
#error Vamp SDK version 2.4 or newer required
#endif

using std::string;

SimpleCepstrum::SimpleCepstrum(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_channels(0),
    m_stepSize(256),
    m_blockSize(1024),
    m_fmin(50),
    m_fmax(1000),
    m_histlen(1),
    m_vflen(1),
    m_clamp(false),
    m_method(InverseSymmetric),
    m_binFrom(0),
    m_binTo(0),
    m_bins(0),
    m_history(0)
{
}

SimpleCepstrum::~SimpleCepstrum()
{
    if (m_history) {
        for (int i = 0; i < m_histlen; ++i) {
            delete[] m_history[i];
        }
        delete[] m_history;
    }
}

string
SimpleCepstrum::getIdentifier() const
{
    return "simple-cepstrum";
}

string
SimpleCepstrum::getName() const
{
    return "Simple Cepstrum";
}

string
SimpleCepstrum::getDescription() const
{
    return "Return simple cepstral data from DFT bins. This plugin is intended for casual inspection of cepstral data. It returns a lot of different sorts of data and is quite slow; it's not a good way to extract a single feature rapidly.";
}

string
SimpleCepstrum::getMaker() const
{
    return "Chris Cannam";
}

int
SimpleCepstrum::getPluginVersion() const
{
    // Increment this each time you release a version that behaves
    // differently from the previous one
    return 1;
}

string
SimpleCepstrum::getCopyright() const
{
    return "Freely redistributable (BSD license)";
}

SimpleCepstrum::InputDomain
SimpleCepstrum::getInputDomain() const
{
    return FrequencyDomain;
}

size_t
SimpleCepstrum::getPreferredBlockSize() const
{
    return 1024;
}

size_t 
SimpleCepstrum::getPreferredStepSize() const
{
    return 256;
}

size_t
SimpleCepstrum::getMinChannelCount() const
{
    return 1;
}

size_t
SimpleCepstrum::getMaxChannelCount() const
{
    return 1;
}

SimpleCepstrum::ParameterList
SimpleCepstrum::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;

    d.identifier = "fmin";
    d.name = "Minimum frequency";
    d.description = "Frequency whose period corresponds to the quefrency of the last cepstrum bin in range";
    d.unit = "Hz";
    d.minValue = m_inputSampleRate / m_blockSize;
    d.maxValue = m_inputSampleRate / 2;
    d.defaultValue = 50;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "fmax";
    d.name = "Maximum frequency";
    d.description = "Frequency whose period corresponds to the quefrency of the first cepstrum bin in range";
    d.unit = "Hz";
    d.minValue = m_inputSampleRate / m_blockSize;
    d.maxValue = m_inputSampleRate / 2;
    d.defaultValue = 1000;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "histlen";
    d.name = "Mean filter history length";
    d.description = "Length of mean filter used for smoothing cepstrum across time bins";
    d.unit = "";
    d.minValue = 1;
    d.maxValue = 10;
    d.defaultValue = 1;
    d.isQuantized = true;
    d.quantizeStep = 1;
    list.push_back(d);

    d.identifier = "vflen";
    d.name = "Vertical filter length";
    d.description = "Length of mean filter used for smoothing cepstrum across quefrency bins";
    d.unit = "";
    d.minValue = 1;
    d.maxValue = 11;
    d.defaultValue = 1;
    d.isQuantized = true;
    d.quantizeStep = 2;
    list.push_back(d);

    d.identifier = "method";
    d.name = "Cepstrum transform method";
    d.unit = "Method to use for calculating cepstrum, starting from the complex short-time Fourier transform of the input audio.\nInverse symmetric - Real part of inverse FFT of log magnitude spectrum, with frequencies above Nyquist reflecting those below it.\nInverse asymmetric - Real part of inverse FFT of log magnitude spectrum, with frequencies above Nyquist set to zero.\nInverse complex - Real part of inverse FFT of complex log spectrum.\nForward magnitude - Magnitude of forward FFT of log magnitude spectrum.\nForward difference - Difference between imaginary and real parts of forward FFT of log magnitude spectrum";
    d.minValue = 0;
    d.maxValue = 4;
    d.defaultValue = 0;
    d.isQuantized = true;
    d.quantizeStep = 1;
    d.valueNames.push_back("Inverse symmetric");
    d.valueNames.push_back("Inverse asymmetric");
    d.valueNames.push_back("Inverse complex");
    d.valueNames.push_back("Forward magnitude");
    d.valueNames.push_back("Forward difference");
    list.push_back(d);

    d.identifier = "clamp";
    d.name = "Clamp negative values in cepstrum at zero";
    d.unit = "If set, no negative values will be returned; they will be replaced by zeroes";
    d.minValue = 0;
    d.maxValue = 1;
    d.defaultValue = 0;
    d.isQuantized = true;
    d.quantizeStep = 1;
    d.valueNames.clear();
    list.push_back(d);

    return list;
}

float
SimpleCepstrum::getParameter(string identifier) const
{
    if (identifier == "fmin") return m_fmin;
    else if (identifier == "fmax") return m_fmax;
    else if (identifier == "histlen") return m_histlen;
    else if (identifier == "vflen") return m_vflen;
    else if (identifier == "clamp") return (m_clamp ? 1 : 0);
    else if (identifier == "method") return (int)m_method;
    else return 0.f;
}

void
SimpleCepstrum::setParameter(string identifier, float value) 
{
    if (identifier == "fmin") m_fmin = value;
    else if (identifier == "fmax") m_fmax = value;
    else if (identifier == "histlen") m_histlen = value;
    else if (identifier == "vflen") m_vflen = value;
    else if (identifier == "clamp") m_clamp = (value > 0.5);
    else if (identifier == "method") m_method = Method(int(value + 0.5));
}

SimpleCepstrum::ProgramList
SimpleCepstrum::getPrograms() const
{
    ProgramList list;
    return list;
}

string
SimpleCepstrum::getCurrentProgram() const
{
    return ""; // no programs
}

void
SimpleCepstrum::selectProgram(string name)
{
}

SimpleCepstrum::OutputList
SimpleCepstrum::getOutputDescriptors() const
{
    OutputList outputs;

    int n = 0;

    OutputDescriptor d;

    d.identifier = "raw_cepstral_peak";
    d.name = "Frequency corresponding to raw cepstral peak";
    d.description = "Return the frequency whose period corresponds to the quefrency with the maximum bin value within the specified range of the cepstrum";
    d.unit = "Hz";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = true;
    d.minValue = m_fmin;
    d.maxValue = m_fmax;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::OneSamplePerStep;
    d.hasDuration = false;
    m_pkOutput = n++;
    outputs.push_back(d);

    d.identifier = "interpolated_peak";
    d.name = "Interpolated peak frequency";
    d.description = "Return the frequency whose period corresponds to the quefrency with the maximum bin value within the specified range of the cepstrum, using parabolic interpolation to estimate the peak quefrency to finer than single bin resolution";
    m_ipkOutput = n++;
    outputs.push_back(d);

    d.identifier = "variance";
    d.name = "Variance of cepstral bins in range";
    d.unit = "";
    d.description = "Return the variance of bin values within the specified range of the cepstrum";
    d.hasKnownExtents = false;
    m_varOutput = n++;
    outputs.push_back(d);

    d.identifier = "peak";
    d.name = "Value at peak";
    d.unit = "";
    d.description = "Return the value found in the maximum-valued bin within the specified range of the cepstrum";
    m_pvOutput = n++;
    outputs.push_back(d);

    d.identifier = "peak_to_rms";
    d.name = "Peak-to-RMS distance";
    d.unit = "";
    d.description = "Return the difference between maximum and root mean square bin values within the specified range of the cepstrum";
    m_p2rOutput = n++;
    outputs.push_back(d);

    d.identifier = "peak_proportion";
    d.name = "Energy around peak";
    d.unit = "";
    d.description = "Return the proportion of total energy that is found in the bins around the peak bin (as far as the nearest local minima), within the specified range of the cepstrum";
    m_ppOutput = n++;
    outputs.push_back(d);

    d.identifier = "peak_to_second_peak";
    d.name = "Peak to second-peak difference";
    d.unit = "";
    d.description = "Return the difference between the value found in the peak bin within the specified range of the cepstrum, and that found in the next highest peak";
    m_pkoOutput = n++;
    outputs.push_back(d);

    d.identifier = "total";
    d.name = "Total energy";
    d.unit = "";
    d.description = "Return the total energy found in all bins within the specified range of the cepstrum";
    m_totOutput = n++;
    outputs.push_back(d);

    d.identifier = "cepstrum";
    d.name = "Cepstrum";
    d.unit = "";
    d.description = "The unprocessed cepstrum bins within the specified range";

    int from = int(m_inputSampleRate / m_fmax);
    int to = int(m_inputSampleRate / m_fmin);
    if (from >= (int)m_blockSize / 2) {
        from = m_blockSize / 2 - 1;
    }
    if (to >= (int)m_blockSize / 2) {
        to = m_blockSize / 2 - 1;
    }
    if (to < from) {
        to = from;
    }
    d.binCount = to - from + 1;
    for (int i = from; i <= to; ++i) {
        float freq = m_inputSampleRate / i;
        char buffer[50];
        sprintf(buffer, "%.2f Hz", freq);
        d.binNames.push_back(buffer);
    }

    d.hasKnownExtents = false;
    m_cepOutput = n++;
    outputs.push_back(d);

    d.identifier = "am";
    d.name = "Cepstrum bins relative to RMS";
    d.description = "The cepstrum bins within the specified range, expressed as a value relative to the root mean square bin value in the range, with values below the RMS clamped to zero";
    m_amOutput = n++;
    outputs.push_back(d);

    d.identifier = "env";
    d.name = "Spectral envelope";
    d.description = "Envelope calculated from the cepstral values below the specified minimum";
    d.binCount = m_blockSize/2 + 1;
    d.binNames.clear();
    for (int i = 0; i < (int)d.binCount; ++i) {
        float freq = (m_inputSampleRate / m_blockSize) * i;
        char buffer[50];
        sprintf(buffer, "%.2f Hz", freq);
        d.binNames.push_back(buffer);
    }
    m_envOutput = n++;
    outputs.push_back(d);

    d.identifier = "es";
    d.name = "Spectrum without envelope";
    d.description = "Magnitude of spectrum values divided by calculated envelope values, to deconvolve the envelope";
    m_esOutput = n++;
    outputs.push_back(d);

    return outputs;
}

bool
SimpleCepstrum::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) return false;

//    std::cerr << "SimpleCepstrum::initialise: channels = " << channels
//	      << ", stepSize = " << stepSize << ", blockSize = " << blockSize
//	      << std::endl;

    m_channels = channels;
    m_stepSize = stepSize;
    m_blockSize = blockSize;

    m_binFrom = int(m_inputSampleRate / m_fmax);
    m_binTo = int(m_inputSampleRate / m_fmin);

    if (m_binFrom >= (int)m_blockSize / 2) {
        m_binFrom = m_blockSize / 2 - 1;
    }
    if (m_binTo >= (int)m_blockSize / 2) {
        m_binTo = m_blockSize / 2 - 1;
    }
    if (m_binTo < m_binFrom) {
        m_binTo = m_binFrom;
    }

    m_bins = (m_binTo - m_binFrom) + 1;

    m_history = new double *[m_histlen];
    for (int i = 0; i < m_histlen; ++i) {
        m_history[i] = new double[m_bins];
    }

    reset();

    return true;
}

void
SimpleCepstrum::reset()
{
    for (int i = 0; i < m_histlen; ++i) {
        for (int j = 0; j < m_bins; ++j) {
            m_history[i][j] = 0.0;
        }
    }
}

void
SimpleCepstrum::filter(const double *cep, double *result)
{
    int hix = m_histlen - 1; // current history index

    // roll back the history
    if (m_histlen > 1) {
        double *oldest = m_history[0];
        for (int i = 1; i < m_histlen; ++i) {
            m_history[i-1] = m_history[i];
        }
        // and stick this back in the newest spot, to recycle
        m_history[hix] = oldest;
    }

    for (int i = 0; i < m_bins; ++i) {
        double v = 0;
        int n = 0;
        // average according to the vertical filter length
        for (int j = -m_vflen/2; j <= m_vflen/2; ++j) {
            int ix = i + m_binFrom + j;
            if (ix >= 0 && ix < (int)m_blockSize) {
                v += cep[ix];
                ++n;
            }
        }
        m_history[hix][i] = v / n;
    }

    for (int i = 0; i < m_bins; ++i) {
        double mean = 0.0;
        for (int j = 0; j < m_histlen; ++j) {
            mean += m_history[j][i];
        }
        mean /= m_histlen;
        result[i] = mean;
    }
}

double
SimpleCepstrum::findInterpolatedPeak(const double *in, int maxbin)
{
    // after jos, 
    // https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html

    if (maxbin < 1 || maxbin > m_bins - 2) {
        return maxbin;
    }

    double alpha = in[maxbin-1];
    double beta  = in[maxbin];
    double gamma = in[maxbin+1];

    double denom = (alpha - 2*beta + gamma);

    if (denom == 0) {
        // flat
        return maxbin;
    }

    double p = ((alpha - gamma) / denom) / 2.0;

    return double(maxbin) + p;
}
   
void
SimpleCepstrum::addStatisticalOutputs(FeatureSet &fs, const double *data)
{
    int n = m_bins;

    double maxval = 0.0;
    int maxbin = 0;

    for (int i = 0; i < n; ++i) {
        if (data[i] > maxval) {
            maxval = data[i];
            maxbin = i;
        }
    }

    double nextPeakVal = 0.0;

    for (int i = 1; i+1 < n; ++i) {
        if (data[i] > data[i-1] &&
            data[i] > data[i+1] &&
            i != maxbin &&
            data[i] > nextPeakVal) {
            nextPeakVal = data[i];
        }
    }

    Feature rf;
    Feature irf;
    if (maxval > 0.0) {
        rf.values.push_back(m_inputSampleRate / (maxbin + m_binFrom));
        double cimax = findInterpolatedPeak(data, maxbin);
        irf.values.push_back(m_inputSampleRate / (cimax + m_binFrom));
    } else {
        rf.values.push_back(0);
        irf.values.push_back(0);
    }
    fs[m_pkOutput].push_back(rf);
    fs[m_ipkOutput].push_back(irf);

    double total = 0;
    for (int i = 0; i < n; ++i) {
        total += data[i];
    }

    Feature tot;
    tot.values.push_back(total);
    fs[m_totOutput].push_back(tot);

    double mean = total / n;

    double totsqr = 0;
    double abstot = 0;
    for (int i = 0; i < n; ++i) {
        totsqr += data[i] * data[i];
        abstot += fabs(data[i]);
    }
    double rms = sqrt(totsqr / n);

    double variance = 0;
    for (int i = 0; i < n; ++i) {
        double dev = fabs(data[i] - mean);
        variance += dev * dev;
    }
    variance /= n;

    double aroundPeak = 0.0;
    double peakProportion = 0.0;
    if (maxval > 0.0) {
        aroundPeak += fabs(maxval);
        int i = maxbin - 1;
        while (i > 0 && data[i] <= data[i+1]) {
            aroundPeak += fabs(data[i]);
            --i;
        }
        i = maxbin + 1;
        while (i < n && data[i] <= data[i-1]) {
            aroundPeak += fabs(data[i]);
            ++i;
        }
    }
    peakProportion = aroundPeak / abstot;
    Feature pp;
    pp.values.push_back(peakProportion);
    fs[m_ppOutput].push_back(pp);

    Feature vf;
    vf.values.push_back(variance);
    fs[m_varOutput].push_back(vf);

    Feature pr;
    pr.values.push_back(maxval - rms);
    fs[m_p2rOutput].push_back(pr);

    Feature pv;
    pv.values.push_back(maxval);
    fs[m_pvOutput].push_back(pv);

    Feature pko;
    if (nextPeakVal != 0.0) {
        pko.values.push_back(maxval - nextPeakVal);
    } else {
        pko.values.push_back(0.0);
    }
    fs[m_pkoOutput].push_back(pko);

    Feature am;
    for (int i = 0; i < n; ++i) {
        if (data[i] < rms) am.values.push_back(0);
        else am.values.push_back(data[i] - rms);
    }
    fs[m_amOutput].push_back(am);
}

void
SimpleCepstrum::addEnvelopeOutputs(FeatureSet &fs, const float *const *inputBuffers, const double *cep)
{
    // Wipe the higher cepstral bins in order to calculate the
    // envelope. This calculation uses the raw cepstrum, not the
    // filtered values (because only values "in frequency range" are
    // filtered).
    int bs = m_blockSize;
    int hs = m_blockSize/2 + 1;

    double *ecep = new double[bs];
    for (int i = 0; i < m_binFrom; ++i) {
        ecep[i] = cep[i] / bs; 
    }
    for (int i = m_binFrom; i < bs; ++i) {
        ecep[i] = 0;
    }
    ecep[0] /= 2;
    if (m_binFrom > 0) {
        ecep[m_binFrom-1] /= 2;
    }

    double *env = new double[bs];
    double *io = new double[bs];

    //!!! This is only right if the previous transform was an inverse one!
    Vamp::FFT::forward(bs, ecep, 0, env, io);

    for (int i = 0; i < hs; ++i) {
        env[i] = exp(env[i]);
    }
    Feature envf;
    for (int i = 0; i < hs; ++i) {
        envf.values.push_back(env[i]);
    }
    fs[m_envOutput].push_back(envf);

    Feature es;
    for (int i = 0; i < hs; ++i) {
        double re = inputBuffers[0][i*2  ] / env[i];
        double im = inputBuffers[0][i*2+1] / env[i];
        double mag = sqrt(re*re + im*im);
        es.values.push_back(mag);
    }
    fs[m_esOutput].push_back(es);

    delete[] env;
    delete[] ecep;
    delete[] io;
}

SimpleCepstrum::FeatureSet
SimpleCepstrum::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    FeatureSet fs;

    int bs = m_blockSize;
    int hs = m_blockSize/2 + 1;

    double *rawcep = new double[bs];
    double *io = new double[bs];

    if (m_method != InverseComplex) {

        double *logmag = new double[bs];
        
        for (int i = 0; i < hs; ++i) {

            double power =
                inputBuffers[0][i*2  ] * inputBuffers[0][i*2  ] +
                inputBuffers[0][i*2+1] * inputBuffers[0][i*2+1];
            double mag = sqrt(power);

            double lm = log(mag + 0.00000001);

            switch (m_method) {
            case InverseSymmetric:
                logmag[i] = lm;
                if (i > 0) logmag[bs - i] = lm;
                break;
            case InverseAsymmetric:
                logmag[i] = lm;
                if (i > 0) logmag[bs - i] = 0;
                break;
            default:
                logmag[bs/2 + i - 1] = lm;
                if (i < hs-1) {
                    logmag[bs/2 - i - 1] = lm;
                }
                break;
            }
        }

        if (m_method == InverseSymmetric ||
            m_method == InverseAsymmetric) {

            Vamp::FFT::inverse(bs, logmag, 0, rawcep, io);

        } else {

            Vamp::FFT::forward(bs, logmag, 0, rawcep, io);

            if (m_method == ForwardDifference) {
                for (int i = 0; i < hs; ++i) {
                    rawcep[i] = fabs(io[i]) - fabs(rawcep[i]);
                }
            } else {
                for (int i = 0; i < hs; ++i) {
                    rawcep[i] = sqrt(rawcep[i]*rawcep[i] + io[i]*io[i]);
                }
            }
        }

        delete[] logmag;

    } else { // InverseComplex

        double *ri = new double[bs];
        double *ii = new double[bs];
        
        for (int i = 0; i < hs; ++i) {
            double re = inputBuffers[0][i*2];
            double im = inputBuffers[0][i*2+1];
            std::complex<double> c(re, im);
            std::complex<double> clog = std::log(c);
            ri[i] = clog.real();
            ii[i] = clog.imag();
            if (i > 0) {
                ri[bs - i] = ri[i];
                ii[bs - i] = -ii[i];
            }
        }

        Vamp::FFT::inverse(bs, ri, ii, rawcep, io);

        delete[] ri;
        delete[] ii;
    }

    if (m_clamp) {
        for (int i = 0; i < bs; ++i) {
            if (rawcep[i] < 0) rawcep[i] = 0;
        }
    }

    delete[] io;

    double *latest = new double[m_bins];
    filter(rawcep, latest);

    int n = m_bins;

    Feature cf;
    for (int i = 0; i < n; ++i) {
        cf.values.push_back(latest[i]);
    }
    fs[m_cepOutput].push_back(cf);

    addStatisticalOutputs(fs, latest);

    addEnvelopeOutputs(fs, inputBuffers, rawcep);

    delete[] latest;
    delete[] rawcep;

    return fs;
}

SimpleCepstrum::FeatureSet
SimpleCepstrum::getRemainingFeatures()
{
    FeatureSet fs;
    return fs;
}
