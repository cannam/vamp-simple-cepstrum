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

#ifndef _SIMPLE_CEPSTRUM_H_
#define _SIMPLE_CEPSTRUM_H_

#include <vamp-sdk/Plugin.h>

class SimpleCepstrum : public Vamp::Plugin
{
public:
    SimpleCepstrum(float inputSampleRate);
    virtual ~SimpleCepstrum();

    std::string getIdentifier() const;
    std::string getName() const;
    std::string getDescription() const;
    std::string getMaker() const;
    int getPluginVersion() const;
    std::string getCopyright() const;

    InputDomain getInputDomain() const;
    size_t getPreferredBlockSize() const;
    size_t getPreferredStepSize() const;
    size_t getMinChannelCount() const;
    size_t getMaxChannelCount() const;

    ParameterList getParameterDescriptors() const;
    float getParameter(std::string identifier) const;
    void setParameter(std::string identifier, float value);

    ProgramList getPrograms() const;
    std::string getCurrentProgram() const;
    void selectProgram(std::string name);

    OutputList getOutputDescriptors() const;

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);
    void reset();

    FeatureSet process(const float *const *inputBuffers,
                       Vamp::RealTime timestamp);

    FeatureSet getRemainingFeatures();

protected:
    size_t m_channels;
    size_t m_stepSize;
    size_t m_blockSize;
    float m_fmin;
    float m_fmax;
    int m_histlen;
    int m_vflen;
    bool m_clamp;

    enum Method {
        InverseSymmetric,
        InverseAsymmetric,
        InverseComplex,
        ForwardMagnitude,
        ForwardDifference
    };

    Method m_method;

    mutable int m_pkOutput;
    mutable int m_ipkOutput;
    mutable int m_varOutput;
    mutable int m_p2rOutput;
    mutable int m_cepOutput;
    mutable int m_pvOutput;
    mutable int m_amOutput;
    mutable int m_envOutput;
    mutable int m_esOutput;
    mutable int m_ppOutput;
    mutable int m_totOutput;
    mutable int m_pkoOutput;

    int m_binFrom;
    int m_binTo;
    int m_bins; // count of "interesting" bins, those returned in m_cepOutput

    double **m_history;
    
    void filter(const double *in, double *out);
    double findInterpolatedPeak(const double *in, int maxbin);
    void addStatisticalOutputs(FeatureSet &fs, const double *data);
    void addEnvelopeOutputs(FeatureSet &fs, const float *const *inputBuffers,
                            const double *raw);
};

#endif
