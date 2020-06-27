#ifndef _PITCH_SHIFTER_WRAP_
#define _PITCH_SHIFTER_WRAP_

#define MAX_FRAME_LENGTH 8192

/*
	���������ο����ͣ�
	http://blogs.zynaptiq.com/bernsee/pitch-shifting-using-the-ft/

	�����ʲ��ޣ�����λ��Ϊ16λ��
*/
class PitchShifter
{
public:
	PitchShifter();
	~PitchShifter();

public:

	/*
		iSemitone �������� -12~12֮�䡣
		|-12		|0|			|12|
		���					���
		iFFFrameSize ����Ҷ�任 �����֡����һ��Ϊ1024,2048��2��n�η�
		iOverSample ��ʱ����Ҷ�任������4�ǱȽϺ��ʵı任������32������á�
		iSampleRate �������Ƶ������
	*/
	void Init(int iSemitone, int iFFTFrameSize, int iOverSample, int iSampleRate);
	
	/*
		psDataIn ��������
		iDataSize psDataIn���ݴ�С
		psDataOut �������

		psDataIn��psDataOut����Ϊͬһ�ڴ档
	*/
	int DoPitchShifter(short* psDataIn, int iDataSize, short* psDataOut);

	void Uninit();

private:
	bool	m_bInit;

	float	m_fInFIFO[MAX_FRAME_LENGTH];
	float	m_fOutFIFO[MAX_FRAME_LENGTH];
	float	m_fFFTworksp[2 * MAX_FRAME_LENGTH];
	float	m_fLastPhase[MAX_FRAME_LENGTH / 2 + 1];
	float	m_fSumPhase[MAX_FRAME_LENGTH / 2 + 1];
	float	m_fOutputAccum[2 * MAX_FRAME_LENGTH];
	float	m_fAnaFreq[MAX_FRAME_LENGTH];
	float	m_fAnaMagn[MAX_FRAME_LENGTH];
	float	m_fSynFreq[MAX_FRAME_LENGTH];
	float	m_fSynMagn[MAX_FRAME_LENGTH];

	float	m_fPitchShifter;//�������
	int		m_iFFTFrameSize;//tfft����֡��
	int		m_iOverSample;//��ʱ����Ҷ�任����
	int		m_iSampleRate;//������Ƶ�Ĳ�����

	int		m_iRover;
};


#endif