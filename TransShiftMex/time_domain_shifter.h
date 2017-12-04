#pragma once

#include <vector>

namespace audapter {

    typedef double dtype;

    class TimeDomainShifter {
    public:
        typedef std::vector<std::pair<dtype, dtype>> PitchShiftSchedule;

        // Constructor of TimeDomainShifter.
        //
        // Args:
        //   sr: Sampling rate.
        //   frameLen: Frame length.
        //   pitchShiftRatio: Pitch-shifting ratio. Suppose the original pitch
        //     value is f0. After shifting, the pitch becomes f0 * ratio. Must
        //     be a positive number.
        //     Therefore:
        //       1.0 corresponds to no shift,
        //       2.0 corresponds to shifting up 1 octave,
        //       0.5 corresponds to shifting down 1 octave,
        //       etc.
        TimeDomainShifter(const int sr,
                          const int frameLen,
                          const PitchShiftSchedule& pitchShiftSchedule);
        virtual ~TimeDomainShifter();

        void reset();

        // Process a frame signal with this instance of TimeDomainShifter.
        // Detects and records zero-crossings (from negative to non-negative).
        // Keep track of the latest pitch cycle.
        //
        // Args:
        //   f: Bandpass filtered signal segment.
        //   x: Original (i.e., not bandpass filtered) signal segment.
        //   y: Output frame pointer. If == nullptr, no output will be generated.
        //   len: Length of the segment. Must match that of len.
        //     TODO(cais): Redundant?
        void processFrame(const dtype* f,
                          const dtype* x,
                          dtype* y,
                          const int len);

        dtype getLatestShiftedPitchHz() const;
        int getLatestInputPitchCycleBegin() const;  // DEBUG
        dtype getDiscontinuity() const;  // DEBUG

    private:
        TimeDomainShifter() {};

        // Copy a segment to the rotating buffer and increment the
        // rotating-buffer pointer.
        //
        // Args:
        //   x: The segment to copy. Assumes that its length is frameLen.
        void copyToRotBuf(const dtype* x);
        const inline dtype accessRotBuf(const int idx) const;
        void checkPitchShiftSchedule() const;
        void checkFrameLen(const int len) const;

        // Generate a pitch-shifted pitch cycle.
        // Uses the last pitch cycle extracted so far.
        // Updates the scratch write index.
        void genShiftedPitchCycle();

        // Determine if the outupt scratch buffer needs a new pitch cycle.
        bool scratchNeedsPitchCycle();

        // Apply smoothing to a part of the scratch.
        //
        // Args:
        //   begin: beginning index of the part to smooth (inclusive).
        //   end: end index of the part to smooth (exclusive).
        void smoothScratch(int begin, int end);

        int getScratchIndex(int i);

        // Get the pitch shift ratio for the current moment, given the
        // prespecified pitch-shift schedule and the current time.
        dtype getPitchShiftRatio() const;

        // methodId:
        //   0: similar alfnie's crop;
        //   1: linear intepolation for b > 0;
        //   2: same as 1, but with smoothing at the junctions of the
        //      pitch-shifted cycles.
        const int methodId = 2;
        const bool verbose = false;
        const dtype maxPitchShiftRatio = 1.5;

        int totalLen;  // Total length received, in # of signal samples.
        int segCount;  // Counter for the number of segments received.
        dtype lastPitchFilteredSample;
        std::vector<int> zcIndices;  // Zero-crossing indices.

        int sr;
        int frameLen;
        PitchShiftSchedule pitchShiftSchedule;

        int bufLen;  // Length of the rotating and scratch buffers.
        dtype* rotBuf;  // Rotating buffer.
        int rotBufPtr;

        int scratchLen;  // Length of the scratch buffer.
        dtype* scratchBuf;  // Output scratch buffer.
        int scratchReadPtr;
        int scratchWritePtr;

        // The beginning idx of the last pitch cycle, inclusive.
        int lastPitchCycleBegin;  
        // The ending idx of the last pitch cycle, exclusive.
        int lastPitchCycleEnd;

        // Experimental. Use 0.0 to deactivate smoothing of pitch cycle length.
        const dtype pitchCycleSmoothingFactor = 0.0;
        dtype trackedPitchCycle;
           
        dtype latestShiftedPitchHz;  // DEBUG
        int latestInputPitchCycleBegin;  // DEBUG For debug only?
        dtype discontinuity;  // DEBUG
    };

};   // namespace Audapter