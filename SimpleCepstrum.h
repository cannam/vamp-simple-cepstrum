/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

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
    bool m_clamp;

    enum Method {
        InverseSymmetric,
        InverseAsymmetric,
        InverseComplex,
        ForwardMagnitude,
        ForwardDifference
    };

    Method m_method;

//    mutable int m_f0Output;
    mutable int m_pkOutput;
    mutable int m_varOutput;
    mutable int m_p2mOutput;
    mutable int m_p2rOutput;
    mutable int m_cepOutput;
    mutable int m_pvOutput;
    mutable int m_amOutput;
    mutable int m_envOutput;
    mutable int m_esOutput;

    int m_binFrom;
    int m_binTo;
    int m_bins; // count of "interesting" bins, those returned in m_cepOutput

    double **m_history;
    
    void filter(const double *in, double *out);
    void fft(unsigned int n, bool inverse,
             double *ri, double *ii, double *ro, double *io);

    void addStatisticalOutputs(FeatureSet &fs, const double *data);
    void addEnvelopeOutputs(FeatureSet &fs, const float *const *inputBuffers,
                            const double *raw);
};

#endif
