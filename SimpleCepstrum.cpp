/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#include "SimpleCepstrum.h"

#include <vector>
#include <algorithm>

#include <cstdio>
#include <cmath>

using std::string;

SimpleCepstrum::SimpleCepstrum(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_channels(0),
    m_stepSize(256),
    m_blockSize(1024),
    m_fmin(50),
    m_fmax(1000),
    m_clamp(false)
{
}

SimpleCepstrum::~SimpleCepstrum()
{
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
    return "Return simple cepstral data from DFT bins";
}

string
SimpleCepstrum::getMaker() const
{
    // Your name here
    return "";
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
    // This function is not ideally named.  It does not necessarily
    // need to say who made the plugin -- getMaker does that -- but it
    // should indicate the terms under which it is distributed.  For
    // example, "Copyright (year). All Rights Reserved", or "GPL"
    return "";
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
    d.description = "";
    d.unit = "Hz";
    d.minValue = m_inputSampleRate / m_blockSize;
    d.maxValue = m_inputSampleRate / 2;
    d.defaultValue = 50;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "fmax";
    d.name = "Maximum frequency";
    d.description = "";
    d.unit = "Hz";
    d.minValue = m_inputSampleRate / m_blockSize;
    d.maxValue = m_inputSampleRate / 2;
    d.defaultValue = 1000;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "clamp";
    d.name = "Clamp negative values in cepstrum at zero";
    d.unit = "";
    d.minValue = 0;
    d.maxValue = 1;
    d.defaultValue = 0;
    d.isQuantized = true;
    d.quantizeStep = 1;
    list.push_back(d);

    return list;
}

float
SimpleCepstrum::getParameter(string identifier) const
{
    if (identifier == "fmin") return m_fmin;
    else if (identifier == "fmax") return m_fmax;
    else if (identifier == "clamp") return (m_clamp ? 1 : 0);
    else return 0.f;
}

void
SimpleCepstrum::setParameter(string identifier, float value) 
{
    if (identifier == "fmin") m_fmin = value;
    else if (identifier == "fmax") m_fmax = value;
    else if (identifier == "clamp") m_clamp = (value > 0.5);
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
    d.identifier = "f0";
    d.name = "Estimated fundamental frequency";
    d.description = "";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = true;
    d.minValue = m_fmin;
    d.maxValue = m_fmax;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::OneSamplePerStep;
    d.hasDuration = false;
    m_f0Output = n++;
    outputs.push_back(d);

    d.identifier = "raw_cepstral_peak";
    d.name = "Frequency corresponding to raw cepstral peak";
    d.unit = "Hz";
    m_rawOutput = n++;
    outputs.push_back(d);

    d.identifier = "variance";
    d.name = "Variance of cepstral bins in range";
    d.unit = "";
    m_varOutput = n++;
    outputs.push_back(d);

    d.identifier = "peak";
    d.name = "Peak value";
    d.unit = "";
    m_pvOutput = n++;
    outputs.push_back(d);

    d.identifier = "peak_to_mean";
    d.name = "Peak-to-mean distance";
    d.unit = "";
    m_p2mOutput = n++;
    outputs.push_back(d);

    d.identifier = "cepstrum";
    d.name = "Cepstrum";
    d.unit = "";

    int from = int(m_inputSampleRate / m_fmax);
    int to = int(m_inputSampleRate / m_fmin);
    if (to >= (int)m_blockSize / 2) {
        to = m_blockSize / 2 - 1;
    }
    d.binCount = to - from + 1;
    for (int i = from; i <= to; ++i) {
        float freq = m_inputSampleRate / i;
        char buffer[10];
        sprintf(buffer, "%.2f", freq);
        d.binNames.push_back(buffer);
    }

    d.hasKnownExtents = false;
    m_cepOutput = n++;
    outputs.push_back(d);

    d.identifier = "am";
    d.name = "Cepstrum bins relative to mean";
    m_amOutput = n++;
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

    return true;
}

void
SimpleCepstrum::reset()
{
}

SimpleCepstrum::FeatureSet
SimpleCepstrum::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    FeatureSet fs;

    int bs = m_blockSize;
    int hs = m_blockSize/2 + 1;

    double *logmag = new double[bs];
    for (int i = 0; i < hs; ++i) {
        double mag = sqrt(inputBuffers[0][i*2  ] * inputBuffers[0][i*2  ] +
                          inputBuffers[0][i*2+1] * inputBuffers[0][i*2+1]);
        logmag[i] = log(mag + 0.000001);
        if (i > 0) {
            logmag[bs - i] = logmag[i];
        }
    }

    double *cep = new double[bs];
    double *discard = new double[bs];
    fft(bs, true, logmag, 0, cep, discard);
    delete[] discard;

    if (m_clamp) {
        for (int i = 0; i < bs; ++i) {
            if (cep[i] < 0) cep[i] = 0;
        }
    }

    int from = int(m_inputSampleRate / m_fmax);
    int to = int(m_inputSampleRate / m_fmin);

    if (to >= bs / 2) {
        to = bs / 2 - 1;
    }

    Feature cf;
    for (int i = from; i <= to; ++i) {
        cf.values.push_back(cep[i]);
    }
    fs[m_cepOutput].push_back(cf);

    float maxval = 0.f;
    int maxbin = 0;

    for (int i = from; i <= to; ++i) {
        if (cep[i] > maxval) {
            maxval = cep[i];
            maxbin = i;
        }
    }

    Feature rf;
    if (maxbin > 0) {
        rf.values.push_back(m_inputSampleRate / maxbin);
    } else {
        rf.values.push_back(0);
    }
    fs[m_rawOutput].push_back(rf);

    float mean = 0;
    for (int i = from; i <= to; ++i) {
        mean += cep[i];
    }
    mean /= (to - from) + 1;

    float variance = 0;
    for (int i = from; i <= to; ++i) {
        float dev = fabsf(cep[i] - mean);
        variance += dev * dev;
    }
    variance /= (to - from) + 1;

    Feature vf;
    vf.values.push_back(variance);
    fs[m_varOutput].push_back(vf);

    Feature pf;
    pf.values.push_back(maxval - mean);
    fs[m_p2mOutput].push_back(pf);

    Feature pv;
    pv.values.push_back(maxval);
    fs[m_pvOutput].push_back(pv);

    Feature am;
    for (int i = from; i <= to; ++i) {
        if (cep[i] < mean) am.values.push_back(0);
        else am.values.push_back(cep[i] - mean);
    }
    fs[m_amOutput].push_back(am);

    delete[] logmag;
    delete[] cep;

    return fs;
}

SimpleCepstrum::FeatureSet
SimpleCepstrum::getRemainingFeatures()
{
    FeatureSet fs;
    return fs;
}

void
SimpleCepstrum::fft(unsigned int n, bool inverse,
                    double *ri, double *ii, double *ro, double *io)
{
    if (!ri || !ro || !io) return;

    unsigned int bits;
    unsigned int i, j, k, m;
    unsigned int blockSize, blockEnd;

    double tr, ti;

    if (n < 2) return;
    if (n & (n-1)) return;

    double angle = 2.0 * M_PI;
    if (inverse) angle = -angle;

    for (i = 0; ; ++i) {
	if (n & (1 << i)) {
	    bits = i;
	    break;
	}
    }

    static unsigned int tableSize = 0;
    static int *table = 0;

    if (tableSize != n) {

	delete[] table;

	table = new int[n];

	for (i = 0; i < n; ++i) {
	
	    m = i;

	    for (j = k = 0; j < bits; ++j) {
		k = (k << 1) | (m & 1);
		m >>= 1;
	    }

	    table[i] = k;
	}

	tableSize = n;
    }

    if (ii) {
	for (i = 0; i < n; ++i) {
	    ro[table[i]] = ri[i];
	    io[table[i]] = ii[i];
	}
    } else {
	for (i = 0; i < n; ++i) {
	    ro[table[i]] = ri[i];
	    io[table[i]] = 0.0;
	}
    }

    blockEnd = 1;

    for (blockSize = 2; blockSize <= n; blockSize <<= 1) {

	double delta = angle / (double)blockSize;
	double sm2 = -sin(-2 * delta);
	double sm1 = -sin(-delta);
	double cm2 = cos(-2 * delta);
	double cm1 = cos(-delta);
	double w = 2 * cm1;
	double ar[3], ai[3];

	for (i = 0; i < n; i += blockSize) {

	    ar[2] = cm2;
	    ar[1] = cm1;

	    ai[2] = sm2;
	    ai[1] = sm1;

	    for (j = i, m = 0; m < blockEnd; j++, m++) {

		ar[0] = w * ar[1] - ar[2];
		ar[2] = ar[1];
		ar[1] = ar[0];

		ai[0] = w * ai[1] - ai[2];
		ai[2] = ai[1];
		ai[1] = ai[0];

		k = j + blockEnd;
		tr = ar[0] * ro[k] - ai[0] * io[k];
		ti = ar[0] * io[k] + ai[0] * ro[k];

		ro[k] = ro[j] - tr;
		io[k] = io[j] - ti;

		ro[j] += tr;
		io[j] += ti;
	    }
	}

	blockEnd = blockSize;
    }
}


