#include "PitchShifter.h"

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <cmath>

#define M_PI 3.14159265358979323846f
#define M_1_PI     0.318309886183790671538f


#define FAST_MATH_TABLE_SIZE  512

const float sinTable_f32[FAST_MATH_TABLE_SIZE + 1] = {
	0.00000000f, 0.01227154f, 0.02454123f, 0.03680722f, 0.04906767f, 0.06132074f,
	0.07356456f, 0.08579731f, 0.09801714f, 0.11022221f, 0.12241068f, 0.13458071f,
	0.14673047f, 0.15885814f, 0.17096189f, 0.18303989f, 0.19509032f, 0.20711138f,
	0.21910124f, 0.23105811f, 0.24298018f, 0.25486566f, 0.26671276f, 0.27851969f,
	0.29028468f, 0.30200595f, 0.31368174f, 0.32531029f, 0.33688985f, 0.34841868f,
	0.35989504f, 0.37131719f, 0.38268343f, 0.39399204f, 0.40524131f, 0.41642956f,
	0.42755509f, 0.43861624f, 0.44961133f, 0.46053871f, 0.47139674f, 0.48218377f,
	0.49289819f, 0.50353838f, 0.51410274f, 0.52458968f, 0.53499762f, 0.54532499f,
	0.55557023f, 0.56573181f, 0.57580819f, 0.58579786f, 0.59569930f, 0.60551104f,
	0.61523159f, 0.62485949f, 0.63439328f, 0.64383154f, 0.65317284f, 0.66241578f,
	0.67155895f, 0.68060100f, 0.68954054f, 0.69837625f, 0.70710678f, 0.71573083f,
	0.72424708f, 0.73265427f, 0.74095113f, 0.74913639f, 0.75720885f, 0.76516727f,
	0.77301045f, 0.78073723f, 0.78834643f, 0.79583690f, 0.80320753f, 0.81045720f,
	0.81758481f, 0.82458930f, 0.83146961f, 0.83822471f, 0.84485357f, 0.85135519f,
	0.85772861f, 0.86397286f, 0.87008699f, 0.87607009f, 0.88192126f, 0.88763962f,
	0.89322430f, 0.89867447f, 0.90398929f, 0.90916798f, 0.91420976f, 0.91911385f,
	0.92387953f, 0.92850608f, 0.93299280f, 0.93733901f, 0.94154407f, 0.94560733f,
	0.94952818f, 0.95330604f, 0.95694034f, 0.96043052f, 0.96377607f, 0.96697647f,
	0.97003125f, 0.97293995f, 0.97570213f, 0.97831737f, 0.98078528f, 0.98310549f,
	0.98527764f, 0.98730142f, 0.98917651f, 0.99090264f, 0.99247953f, 0.99390697f,
	0.99518473f, 0.99631261f, 0.99729046f, 0.99811811f, 0.99879546f, 0.99932238f,
	0.99969882f, 0.99992470f, 1.00000000f, 0.99992470f, 0.99969882f, 0.99932238f,
	0.99879546f, 0.99811811f, 0.99729046f, 0.99631261f, 0.99518473f, 0.99390697f,
	0.99247953f, 0.99090264f, 0.98917651f, 0.98730142f, 0.98527764f, 0.98310549f,
	0.98078528f, 0.97831737f, 0.97570213f, 0.97293995f, 0.97003125f, 0.96697647f,
	0.96377607f, 0.96043052f, 0.95694034f, 0.95330604f, 0.94952818f, 0.94560733f,
	0.94154407f, 0.93733901f, 0.93299280f, 0.92850608f, 0.92387953f, 0.91911385f,
	0.91420976f, 0.90916798f, 0.90398929f, 0.89867447f, 0.89322430f, 0.88763962f,
	0.88192126f, 0.87607009f, 0.87008699f, 0.86397286f, 0.85772861f, 0.85135519f,
	0.84485357f, 0.83822471f, 0.83146961f, 0.82458930f, 0.81758481f, 0.81045720f,
	0.80320753f, 0.79583690f, 0.78834643f, 0.78073723f, 0.77301045f, 0.76516727f,
	0.75720885f, 0.74913639f, 0.74095113f, 0.73265427f, 0.72424708f, 0.71573083f,
	0.70710678f, 0.69837625f, 0.68954054f, 0.68060100f, 0.67155895f, 0.66241578f,
	0.65317284f, 0.64383154f, 0.63439328f, 0.62485949f, 0.61523159f, 0.60551104f,
	0.59569930f, 0.58579786f, 0.57580819f, 0.56573181f, 0.55557023f, 0.54532499f,
	0.53499762f, 0.52458968f, 0.51410274f, 0.50353838f, 0.49289819f, 0.48218377f,
	0.47139674f, 0.46053871f, 0.44961133f, 0.43861624f, 0.42755509f, 0.41642956f,
	0.40524131f, 0.39399204f, 0.38268343f, 0.37131719f, 0.35989504f, 0.34841868f,
	0.33688985f, 0.32531029f, 0.31368174f, 0.30200595f, 0.29028468f, 0.27851969f,
	0.26671276f, 0.25486566f, 0.24298018f, 0.23105811f, 0.21910124f, 0.20711138f,
	0.19509032f, 0.18303989f, 0.17096189f, 0.15885814f, 0.14673047f, 0.13458071f,
	0.12241068f, 0.11022221f, 0.09801714f, 0.08579731f, 0.07356456f, 0.06132074f,
	0.04906767f, 0.03680722f, 0.02454123f, 0.01227154f, 0.00000000f, -0.01227154f,
	-0.02454123f, -0.03680722f, -0.04906767f, -0.06132074f, -0.07356456f,
	-0.08579731f, -0.09801714f, -0.11022221f, -0.12241068f, -0.13458071f,
	-0.14673047f, -0.15885814f, -0.17096189f, -0.18303989f, -0.19509032f,
	-0.20711138f, -0.21910124f, -0.23105811f, -0.24298018f, -0.25486566f,
	-0.26671276f, -0.27851969f, -0.29028468f, -0.30200595f, -0.31368174f,
	-0.32531029f, -0.33688985f, -0.34841868f, -0.35989504f, -0.37131719f,
	-0.38268343f, -0.39399204f, -0.40524131f, -0.41642956f, -0.42755509f,
	-0.43861624f, -0.44961133f, -0.46053871f, -0.47139674f, -0.48218377f,
	-0.49289819f, -0.50353838f, -0.51410274f, -0.52458968f, -0.53499762f,
	-0.54532499f, -0.55557023f, -0.56573181f, -0.57580819f, -0.58579786f,
	-0.59569930f, -0.60551104f, -0.61523159f, -0.62485949f, -0.63439328f,
	-0.64383154f, -0.65317284f, -0.66241578f, -0.67155895f, -0.68060100f,
	-0.68954054f, -0.69837625f, -0.70710678f, -0.71573083f, -0.72424708f,
	-0.73265427f, -0.74095113f, -0.74913639f, -0.75720885f, -0.76516727f,
	-0.77301045f, -0.78073723f, -0.78834643f, -0.79583690f, -0.80320753f,
	-0.81045720f, -0.81758481f, -0.82458930f, -0.83146961f, -0.83822471f,
	-0.84485357f, -0.85135519f, -0.85772861f, -0.86397286f, -0.87008699f,
	-0.87607009f, -0.88192126f, -0.88763962f, -0.89322430f, -0.89867447f,
	-0.90398929f, -0.90916798f, -0.91420976f, -0.91911385f, -0.92387953f,
	-0.92850608f, -0.93299280f, -0.93733901f, -0.94154407f, -0.94560733f,
	-0.94952818f, -0.95330604f, -0.95694034f, -0.96043052f, -0.96377607f,
	-0.96697647f, -0.97003125f, -0.97293995f, -0.97570213f, -0.97831737f,
	-0.98078528f, -0.98310549f, -0.98527764f, -0.98730142f, -0.98917651f,
	-0.99090264f, -0.99247953f, -0.99390697f, -0.99518473f, -0.99631261f,
	-0.99729046f, -0.99811811f, -0.99879546f, -0.99932238f, -0.99969882f,
	-0.99992470f, -1.00000000f, -0.99992470f, -0.99969882f, -0.99932238f,
	-0.99879546f, -0.99811811f, -0.99729046f, -0.99631261f, -0.99518473f,
	-0.99390697f, -0.99247953f, -0.99090264f, -0.98917651f, -0.98730142f,
	-0.98527764f, -0.98310549f, -0.98078528f, -0.97831737f, -0.97570213f,
	-0.97293995f, -0.97003125f, -0.96697647f, -0.96377607f, -0.96043052f,
	-0.95694034f, -0.95330604f, -0.94952818f, -0.94560733f, -0.94154407f,
	-0.93733901f, -0.93299280f, -0.92850608f, -0.92387953f, -0.91911385f,
	-0.91420976f, -0.90916798f, -0.90398929f, -0.89867447f, -0.89322430f,
	-0.88763962f, -0.88192126f, -0.87607009f, -0.87008699f, -0.86397286f,
	-0.85772861f, -0.85135519f, -0.84485357f, -0.83822471f, -0.83146961f,
	-0.82458930f, -0.81758481f, -0.81045720f, -0.80320753f, -0.79583690f,
	-0.78834643f, -0.78073723f, -0.77301045f, -0.76516727f, -0.75720885f,
	-0.74913639f, -0.74095113f, -0.73265427f, -0.72424708f, -0.71573083f,
	-0.70710678f, -0.69837625f, -0.68954054f, -0.68060100f, -0.67155895f,
	-0.66241578f, -0.65317284f, -0.64383154f, -0.63439328f, -0.62485949f,
	-0.61523159f, -0.60551104f, -0.59569930f, -0.58579786f, -0.57580819f,
	-0.56573181f, -0.55557023f, -0.54532499f, -0.53499762f, -0.52458968f,
	-0.51410274f, -0.50353838f, -0.49289819f, -0.48218377f, -0.47139674f,
	-0.46053871f, -0.44961133f, -0.43861624f, -0.42755509f, -0.41642956f,
	-0.40524131f, -0.39399204f, -0.38268343f, -0.37131719f, -0.35989504f,
	-0.34841868f, -0.33688985f, -0.32531029f, -0.31368174f, -0.30200595f,
	-0.29028468f, -0.27851969f, -0.26671276f, -0.25486566f, -0.24298018f,
	-0.23105811f, -0.21910124f, -0.20711138f, -0.19509032f, -0.18303989f,
	-0.17096189f, -0.15885814f, -0.14673047f, -0.13458071f, -0.12241068f,
	-0.11022221f, -0.09801714f, -0.08579731f, -0.07356456f, -0.06132074f,
	-0.04906767f, -0.03680722f, -0.02454123f, -0.01227154f, -0.00000000f
};

float fastsin(float x)
{
	float sinVal, fract, in;                           /* Temporary variables for input, output */
	unsigned short  index;                                        /* Index variable */
	float a, b;                                        /* Two nearest output values */
	int n;
	float findex;

	/* input x is in radians */
	/* Scale the input to [0 1] range from [0 2*PI] , divide input by 2*pi */
	in = x * 0.159154943092f;

	/* Calculation of floor value of input */
	n = (int)in;

	/* Make negative values towards -infinity */
	if (x < 0.0f)
	{
		n--;
	}

	/* Map input value to [0 1] */
	in = in - (float)n;

	/* Calculation of index of the table */
	findex = (float)FAST_MATH_TABLE_SIZE * in;
	if (findex >= 512.0f) {
		findex -= 512.0f;
	}

	index = ((unsigned short)findex) & 0x1ff;

	/* fractional value calculation */
	fract = findex - (float)index;

	/* Read two nearest values of input value from the sin table */
	a = sinTable_f32[index];
	b = sinTable_f32[index + 1];

	/* Linear interpolation process */
	sinVal = (1.0f - fract)*a + fract*b;

	/* Return the output value */
	return (sinVal);
}

float  fastcos(float x)
{
	float cosVal, fract, in;                   /* Temporary variables for input, output */
	unsigned short index;                                /* Index variable */
	float a, b;                                /* Two nearest output values */
	int n;
	float findex;

	/* input x is in radians */
	/* Scale the input to [0 1] range from [0 2*PI] , divide input by 2*pi, add 0.25 (pi/2) to read sine table */
	in = x * 0.159154943092f + 0.25f;

	/* Calculation of floor value of input */
	n = (int)in;

	/* Make negative values towards -infinity */
	if (in < 0.0f)
	{
		n--;
	}

	/* Map input value to [0 1] */
	in = in - (float)n;

	/* Calculation of index of the table */
	findex = (float)FAST_MATH_TABLE_SIZE * in;
	index = ((unsigned short)findex) & 0x1ff;

	/* fractional value calculation */
	fract = findex - (float)index;

	/* Read two nearest values of input value from the cos table */
	a = sinTable_f32[index];
	b = sinTable_f32[index + 1];

	/* Linear interpolation process */
	cosVal = (1.0f - fract)*a + fract*b;

	/* Return the output value */
	return (cosVal);
}

float fastAtan2(float y, float x) {
	//float coeff_1 = PI/4.0f;
	float coeff_1 = 0.785398163f;
	float coeff_2 = 3.0f * coeff_1;
	float abs_y = abs(y) + 0.00000001f;    // kludge to prevent 0/0 condition;
	float angle;
	if (x >= 0.0f) {
		float r = (x - abs_y) / (x + abs_y);
		angle = coeff_1 - coeff_1 * r;
	}
	else {
		float r = (x + abs_y) / (abs_y - x);
		angle = coeff_2 - coeff_1 * r;
	}
	return y < 0.0f ? -angle : angle;
}

float fastSqrt(float v)
{
	int i;
	float x2, y;
	const float threehalfs = 1.5F;

	x2 = v * 0.5F;
	y = v;
	i = *(int *)&y;
	i = 0x5f375a86 - (i >> 1);
	y = *(float *)&i;
	y = y * (threehalfs - (x2 * y * y));
	y = y * (threehalfs - (x2 * y * y));
	y = y * (threehalfs - (x2 * y * y));
	return v*y;
}

/*

12/12/02, smb

PLEASE NOTE:

There have been some reports on domain errors when the atan2() function was used
as in the above code. Usually, a domain error should not interrupt the program flow
(maybe except in Debug mode) but rather be handled "silently" and a global variable
should be set according to this error. However, on some occasions people ran into
this kind of scenario, so a replacement atan2() function is provided here.

If you are experiencing domain errors and your program stops, simply replace all
instances of atan2() with calls to the smbAtan2() function below.

*/
float smbAtan2(float x, float y)
{
	float signx;
	if (x > 0.0f) signx = 1.0f;
	else signx = -1.0f;

	if (x == 0.0f) return 0.0f;
	if (y == 0.0f) return signx * M_PI / 2.0f;

	return fastAtan2(x, y);
}
// -----------------------------------------------------------------------------------------------------------------

/*
FFT routine, (C)1996 S.M.Bernsee. Sign = -1 is FFT, 1 is iFFT (inverse)
Fills fftBuffer[0...2*fftFrameSize-1] with the Fourier transform of the
time domain data in fftBuffer[0...2*fftFrameSize-1]. The FFT array takes
and returns the cosine and sine parts in an interleaved manner, ie.
fftBuffer[0] = cosPart[0], fftBuffer[1] = sinPart[0], asf. fftFrameSize
must be a power of 2. It expects a complex input signal (see footnote 2),
ie. when working with 'common' audio signals our input signal has to be
passed as {in[0],0.,in[1],0.,in[2],0.,...} asf. In that case, the transform
of the frequencies of interest is in fftBuffer[0...fftFrameSize].
*/
void smbFft(float *fftBuffer, long fftFrameSize, long sign)
{
	float wr, wi, arg, *p1, *p2, temp;
	float tr, ti, ur, ui, *p1r, *p1i, *p2r, *p2i;
	long i, bitm, j, le, le2, k;

	for (i = 2; i < 2 * fftFrameSize - 2; i += 2) {
		for (bitm = 2, j = 0; bitm < 2 * fftFrameSize; bitm <<= 1) {
			if (i & bitm) j++;
			j <<= 1;
		}
		if (i < j) {
			p1 = fftBuffer + i; p2 = fftBuffer + j;
			temp = *p1; *(p1++) = *p2;
			*(p2++) = temp; temp = *p1;
			*p1 = *p2; *p2 = temp;
		}
	}
	long k_size = (long)(log((float)fftFrameSize) / log(2.0f) + 0.5f);
	for (k = 0, le = 2; k < k_size; k++) {

		le <<= 1;
		le2 = le >> 1;
		ur = 1.0f;
		ui = 0.0f;
		arg = M_PI / (le2 >> 1);
		wr = fastcos(arg);
		wi = sign*fastsin(arg);
		for (j = 0; j < le2; j += 2) {
			p1r = fftBuffer + j; p1i = p1r + 1;
			p2r = p1r + le2; p2i = p2r + 1;
			for (i = j; i < 2 * fftFrameSize; i += le) {
				tr = *p2r * ur - *p2i * ui;
				ti = *p2r * ui + *p2i * ur;
				*p2r = *p1r - tr; *p2i = *p1i - ti;
				*p1r += tr; *p1i += ti;
				p1r += le; p1i += le;
				p2r += le; p2i += le;
			}
			tr = ur*wr - ui*wi;
			ui = ur*wi + ui*wr;
			ur = tr;
		}
	}
}

PitchShifter::PitchShifter()
	: m_bInit(false)
	, m_iRover(0)
	, m_fPitchShifter(1.0)
	, m_iFFTFrameSize(1024)
	, m_iOverSample(4)
	, m_iSampleRate(48000)
{

}

PitchShifter::~PitchShifter()
{

}


void PitchShifter::Init(int iSemitone, int iFFTFrameSize, int iOverSample, int sampleRate)
{
	if (!m_bInit)
	{
		memset(m_fInFIFO, 0, MAX_FRAME_LENGTH * sizeof(float));
		memset(m_fOutFIFO, 0, MAX_FRAME_LENGTH * sizeof(float));
		memset(m_fFFTworksp, 0, 2 * MAX_FRAME_LENGTH * sizeof(float));
		memset(m_fLastPhase, 0, (MAX_FRAME_LENGTH / 2 + 1) * sizeof(float));
		memset(m_fSumPhase, 0, (MAX_FRAME_LENGTH / 2 + 1) * sizeof(float));
		memset(m_fOutputAccum, 0, 2 * MAX_FRAME_LENGTH * sizeof(float));
		memset(m_fAnaFreq, 0, MAX_FRAME_LENGTH * sizeof(float));
		memset(m_fAnaMagn, 0, MAX_FRAME_LENGTH * sizeof(float));

		m_iRover = 0;
		m_fPitchShifter = pow(2.0f, iSemitone / 12.0f);
		m_iOverSample = iOverSample;
		m_iFFTFrameSize = iFFTFrameSize;
		m_iSampleRate = sampleRate;

		printf("m_fPitchShifter:%g.\n", m_fPitchShifter);
		m_bInit = true;
	}
}


int PitchShifter::DoPitchShifter(short* psDataIn, int iDataSize, short* psDataOut)
{
	float fPhase, tmp, real, imag;

	long qpd, index;
	float M_2PI = 2.0f * M_PI;
	float W_2PI = M_2PI / m_iFFTFrameSize;
	/* set up some handy variables */
	long fftFrameSize2 = m_iFFTFrameSize / 2;
	long stepSize = m_iFFTFrameSize / m_iOverSample;
	float freqPerBin = m_iSampleRate / (float)m_iFFTFrameSize;
	float expct = W_2PI*(float)stepSize;
	float pi_osamp = (M_2PI / m_iOverSample);
	float pi_osamp_expct = pi_osamp + expct;
	float pi_osamp_freqPerBin = pi_osamp / freqPerBin;
	float osamp_freqPerBin = m_iOverSample * freqPerBin / M_2PI;
	float fft_osamp = 1.0f / (fftFrameSize2 * m_iOverSample);

	long inFifoLatency = m_iFFTFrameSize - stepSize;
	if (m_iRover == 0)
		m_iRover = inFifoLatency;

	/* main processing loop */
	int iSampleSum = iDataSize / 2;//总的采样数
	for (int i = 0; i < iSampleSum; i++)
	{
		/* As long as we have not yet collected enough data just read in */
		m_fInFIFO[m_iRover] = 1.0 * psDataIn[i] / 32768;
		psDataOut[i] = m_fOutFIFO[m_iRover - inFifoLatency] * 32768;
		m_iRover++;

		/* now we have enough data for processing */
		if (m_iRover >= m_iFFTFrameSize) {
			m_iRover = inFifoLatency;

			/* do windowing and re,im interleave */
			float * sFFTworksp = m_fFFTworksp;
			for (int k = 0; k < m_iFFTFrameSize; k++) {
				sFFTworksp[0] = m_fInFIFO[k] * (-0.5f*fastcos(W_2PI*k) + 0.5f);
				sFFTworksp[1] = 0.0f;
				sFFTworksp += 2;
			}
			/* ***************** ANALYSIS ******************* */
			/* do transform */
			smbFft(m_fFFTworksp, m_iFFTFrameSize, -1);

			/* this is the analysis step */
			float *pFFTworksp = m_fFFTworksp;
			for (int k = 0; k <= fftFrameSize2; k++) {
				/* de-interlace FFT buffer */
				real = pFFTworksp[0];
				imag = pFFTworksp[1];
				/* compute magnitude and phase */
				m_fAnaMagn[k] = 2.0f * fastSqrt(real*real + imag*imag);
				fPhase = fastAtan2(imag, real);
				/* compute phase difference */
				/* subtract expected phase difference */
				tmp = (fPhase - m_fLastPhase[k]) - (float)k*expct;
				m_fLastPhase[k] = fPhase;

				/* map delta phase into +/- Pi interval */
				qpd = tmp * M_1_PI;
				if (qpd >= 0)
					qpd += qpd & 1;
				else
					qpd -= qpd & 1;
				/* get deviation from bin frequency from the +/- Pi interval */
				/* compute the k-th partials' true frequency */
				/* store magnitude and true frequency in analysis arrays */
				m_fAnaFreq[k] = (k + (tmp - M_PI * qpd)  * osamp_freqPerBin);
				pFFTworksp += 2;
			}

			/* ***************** PROCESSING ******************* */
			/* this does the actual pitch shifting */
			memset(m_fSynMagn, 0, m_iFFTFrameSize * sizeof(float));
			memset(m_fSynFreq, 0, m_iFFTFrameSize * sizeof(float));
			for (int k = 0; k <= fftFrameSize2; k++) {
				index = k * m_fPitchShifter;
				if (index <= fftFrameSize2) {
					m_fSynMagn[index] += m_fAnaMagn[k];
					m_fSynFreq[index] = m_fAnaFreq[k] * m_fPitchShifter;
				}
			}

			/* ***************** SYNTHESIS ******************* */
			/* this is the synthesis step */
			float *tFFTworksp = m_fFFTworksp;
			for (int k = 0; k <= fftFrameSize2; k++) {
				m_fSumPhase[k] += pi_osamp_freqPerBin * m_fSynFreq[k] - pi_osamp_expct*k;
				tFFTworksp[0] = m_fSynMagn[k] * fastcos(m_fSumPhase[k]);
				tFFTworksp[1] = m_fSynMagn[k] * fastsin(m_fSumPhase[k]);
				tFFTworksp += 2;
			}

			/* zero negative frequencies */
			memset(m_fFFTworksp + (m_iFFTFrameSize + 2), 0, sizeof(float) * 2 * m_iFFTFrameSize);
			/* do inverse transform */
			smbFft(m_fFFTworksp, m_iFFTFrameSize, 1);
			/* do windowing and add to output accumulator */
			float *ppFFTworksp = m_fFFTworksp;
			for (int k = 0; k < m_iFFTFrameSize; k++) {
				m_fOutputAccum[k] += (-fastcos(W_2PI * k) *ppFFTworksp[0] * fft_osamp) + ppFFTworksp[0] * fft_osamp;
				ppFFTworksp += 2;
			}
			memcpy(m_fOutFIFO, m_fOutputAccum, stepSize * sizeof(float));
			/* shift accumulator */
			memmove(m_fOutputAccum, m_fOutputAccum + stepSize, m_iFFTFrameSize * sizeof(float));
			/* move input FIFO */
			memcpy(m_fInFIFO, m_fInFIFO + stepSize, sizeof(float)*inFifoLatency);
		}
	}

	return 0;
}

void PitchShifter::Uninit()
{
	//nothing to do
	m_bInit = false;
}
