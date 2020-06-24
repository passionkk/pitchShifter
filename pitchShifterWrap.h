#ifndef _PITCH_SHIFTER_WRAP_
#define _PITCH_SHIFTER_WRAP_


class PitchShifter
{
public:
	PitchShifter();
	~PitchShifter();

public:
	void Init();
	/*
		fPitchShift ������Ҫ����
		iSampleSum �ܵĲ�����
		iFFFrameSize ���ٸ���Ҷ �����֡����һ��Ϊ1024,2048��2��n�η�
		oversampling factor
	*/
	int DoPitchShifter(float fPitchShift, int iSampleSum, int iFFTFrameSize, long osamp, float sampleRate, short* indata, short *outdata);
	void Uninit();
};


#endif