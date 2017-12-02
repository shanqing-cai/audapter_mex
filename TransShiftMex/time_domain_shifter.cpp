#include <sstream>
#include <math.h>
#include <vector>
#include "mex.h"

#include "DSPF.h"
#include "time_domain_shifter.h"
#include "utils.h"

namespace audapter {

    TimeDomainShifter::TimeDomainShifter(
        const int sr,
        const int frameLen,
        const PitchShiftSchedule& pitchShiftSchedule) :
        sr(sr), frameLen(frameLen), pitchShiftSchedule(pitchShiftSchedule),
        rotBufPtr(0), scratchReadPtr(0), scratchWritePtr(0),
        latestShiftedPitchHz(0.0) {
        checkPitchShiftSchedule();

        bufLen = static_cast<int>(static_cast<dtype>(sr) / 25.0);
        bufLen = frameLen * (static_cast<int>(bufLen / frameLen) + 1);

        rotBuf = new dtype[bufLen];
        scratchLen = 2 * bufLen;
        scratchBuf = new dtype[scratchLen];
        reset();

        if (verbose) {
            mexPrintf("bufLen = %d; scratchLen = %d\n", bufLen, scratchLen);
        }
    }

    TimeDomainShifter::~TimeDomainShifter() {
        if (rotBuf) {
            delete[] rotBuf;
        }
        if (scratchBuf) {
            delete[] scratchBuf;
        }
    }

    void TimeDomainShifter::reset() {
        totalLen = 0;
        segCount = 0;        
        lastPitchFilteredSample = 0.0;
        zcIndices.clear();

        for (int i = 0; i < bufLen; ++i) {
            rotBuf[i] = 0.0;
        }
        rotBufPtr = 0;

        for (int i = 0; i < scratchLen; ++i) {
            scratchBuf[i] = 0.0;
        }
        scratchReadPtr = 0;
        scratchWritePtr = 0;

        lastPitchCycleBegin = 0;
        lastPitchCycleEnd = 0;
    }

    dtype TimeDomainShifter::getLatestShiftedPitchHz() const {
        return latestShiftedPitchHz;
    }

    void TimeDomainShifter::processFrame(const dtype* f,
                                         const dtype* x,
                                         dtype* y,
                                         const int len) {
        checkFrameLen(len);
        const int zcCountOld = static_cast<int>(zcIndices.size());

        if (segCount > 0) {
            if (lastPitchFilteredSample < 0.0 && f[0] >= 0.0) {
                zcIndices.push_back(totalLen);
            }
        }
        for (int i = 1; i < len; ++i) {
            if (f[i - 1] < 0.0 && f[i] > 0.0) {
                zcIndices.push_back(totalLen + i);
            }
        }

        // Write the segment to the rotating buffer.
        copyToRotBuf(x);

        // Extract last pitch cycle.
        if (zcIndices.size() > zcCountOld) {
            // There is a new pitch cycle. Extract it.
            lastPitchCycleBegin =
                zcIndices.size() > 1 ? zcIndices[zcIndices.size() - 2] : 0;
            lastPitchCycleEnd = zcIndices[zcIndices.size() - 1];

            if (verbose) {
                mexPrintf("ps0 =[");
                for (int i = lastPitchCycleBegin; i < lastPitchCycleEnd; ++i) {
                    mexPrintf("%f", accessRotBuf(i));
                    if (i < lastPitchCycleEnd - 1) {
                        mexPrintf(",");
                    }
                }
                mexPrintf("]\n");
            }
        }

        while (scratchNeedsPitchCycle()) {
            // Use a while loop, because it is possible that more than one
            // additional pitch cycle will be required to fill the next
            // output frame.
            genShiftedPitchCycle();
        }

        if (y) {  // Generate output.
            if (zcIndices.empty()) {
                // No zero-crossing has been detected yet, which means no pitch
                // cycle has been extracted yet.
                // DSPF_dp_blk_move(x, y, len);
                for (int i = 0; i < len; ++i) {
                    y[0] = 0.0;
                }
            } else {
                for (int i = 0; i < frameLen; ++i) {
                    y[i] = scratchBuf[(scratchReadPtr + i) % scratchLen];
                }
                scratchReadPtr = (scratchReadPtr + frameLen) % scratchLen;
            }
        }

        segCount++;
        lastPitchFilteredSample = f[len - 1];
        totalLen += len;
    }

    void TimeDomainShifter::copyToRotBuf(const dtype* x) {
        DSPF_dp_blk_move(x, rotBuf + rotBufPtr, frameLen);
        rotBufPtr = (rotBufPtr + frameLen) % bufLen;
    }

    const inline dtype TimeDomainShifter::accessRotBuf(const int idx) const {
        const int nonNegIdx = idx >= 0 ? idx : idx + bufLen;
        return rotBuf[nonNegIdx % bufLen];
    }

    // TODO(cais): Remove.
    /*const inline dtype TimeDomainShifter::accessScratchBuf(
        const int idx) const {
        const int nonNegIdx = idx >= 0 ? idx : idx + scratchLen;
        return scratchBuf[nonNegIdx % bufLen];
    }*/

    bool TimeDomainShifter::scratchNeedsPitchCycle() {
        if (zcIndices.empty()) {
            return false;
        }

        int unwrappedWritePtr = scratchWritePtr;
        if (scratchReadPtr > scratchWritePtr) {
            // There is wrapping around.
            unwrappedWritePtr += scratchLen;
        }
        return scratchReadPtr + frameLen >= unwrappedWritePtr;
    }

    dtype TimeDomainShifter::getPitchShiftRatio() const {
        if (pitchShiftSchedule.empty()) {
            // The pitch-shift schedule is empty: no shifting.
            return 1.0;
        }

        const dtype nowSec = static_cast<dtype>(totalLen) / sr;
        // Find the last interval that nowSec falls into.
        if (nowSec >= pitchShiftSchedule[pitchShiftSchedule.size() - 1].first) {
            return pitchShiftSchedule[pitchShiftSchedule.size() - 1].second;
        }
        else {
            for (size_t i = pitchShiftSchedule.size() - 1; i > 0; --i) {
                if (nowSec >= pitchShiftSchedule[i - 1].first &&
                    nowSec < pitchShiftSchedule[i].first) {
                    // Do linear interpolation.
                    const dtype timeFrac =
                        (nowSec - pitchShiftSchedule[i - 1].first) /
                        (pitchShiftSchedule[i].first - pitchShiftSchedule[i - 1].first);
                    return pitchShiftSchedule[i - 1].second + timeFrac * (
                        pitchShiftSchedule[i].second - pitchShiftSchedule[i - 1].second);
                }
            }
        }
    }

    void TimeDomainShifter::genShiftedPitchCycle() {
        // Length of the latest unshifted pitch cycle.
        const int n0 = lastPitchCycleEnd - lastPitchCycleBegin;
        const dtype pitchShiftRatio = getPitchShiftRatio();
        latestShiftedPitchHz = static_cast<dtype>(sr) / n0 * pitchShiftRatio;   
        // Length of the shifted pitch cycle.
        //mexPrintf("pitchShiftRatio = %f\n", pitchShiftRatio);  // DEBUG
        const int n1 = static_cast<int>(
            round(static_cast<dtype>(n0) / pitchShiftRatio));

        // First, zero out the n1 samples in scratch buffer.
        for (int i = 0; i < n1; ++i) {
            scratchBuf[(scratchWritePtr + i) % scratchLen] = 0.0;
        }

        int b = 0;
        dtype offset = 0.0;
        while (b < n0) {
            if (verbose) {
                mexPrintf("ps_b = [");
            }
            for (int i = 0; i < n1; ++i) {
                dtype added = 0.0;
                int accessIdx = lastPitchCycleBegin + b + i;
                if (accessIdx > lastPitchCycleEnd - 1) {
                    accessIdx = lastPitchCycleEnd - 1;
                }
                added = accessRotBuf(accessIdx);
                if (methodId != 0) {
                    if (b > 0) {
                        const dtype y2 = accessRotBuf(lastPitchCycleBegin + b);
                        const dtype y3 = accessRotBuf(lastPitchCycleEnd - 1);
                        const dtype y2prime = 0.5 * (y2 - y3);
                        const dtype y3prime = 0.5 * (y3 - y2);
                        added = y2prime + (y3prime - y2prime) * static_cast<dtype>(i)
                            / static_cast<dtype>(n1 - 1);
                    }
                }
                scratchBuf[(scratchWritePtr + i) % scratchLen] += added;
                if (methodId == 0) {
                    offset += added;
                }

                if (verbose) {
                    mexPrintf("%f", added);
                    if (i < n1 - 1) {
                        mexPrintf(",");
                    }
                }
            }
            if (verbose) {
                mexPrintf("]\n");
            }
            b += n1;
        }

        if (methodId == 0) {
            // TODO(cais): Optimize: The baseline can be calculated beforehand.
            dtype baseline = 0.0;
            for (int i = lastPitchCycleBegin; i < lastPitchCycleEnd; ++i) {
                baseline += accessRotBuf(i);
            }
            offset -= baseline;
            offset /= (lastPitchCycleEnd - lastPitchCycleBegin);

            // Substract the offset.
            for (int i = 0; i < n1; ++i) {
                scratchBuf[(scratchWritePtr + i) % scratchLen] -= offset;
            }
        }

        if (methodId == 2) {
            smoothScratch(scratchWritePtr - 1, scratchWritePtr + 2);
        }

        if (verbose) {
            mexPrintf("offset = %f\n", offset);
            mexPrintf("ps_c = [");
            for (int i = 0; i < n1; ++i) {
                mexPrintf("%f", scratchBuf[(scratchWritePtr + i) % scratchLen]);
                if (i < n1 - 1) {
                    mexPrintf(",");
                }
            }
            mexPrintf("]\n");
        }

        // Update scratchWritePtr.
        scratchWritePtr += n1;
        if (scratchWritePtr >= scratchLen) {
            scratchWritePtr %= scratchLen;
        }
    }

    void TimeDomainShifter::checkPitchShiftSchedule() const {
        // If pitchShfitSchedule is empty, return right away. (It means no pitch
        // shifting.)
        if (pitchShiftSchedule.empty()) {
            return;
        }

        // Check that the first time point is 0.0.
        if (pitchShiftSchedule[0].first != 0.0) {
            std::ostringstream errMsg;
            errMsg << "The first time point in the pitch-shift schedule "
                << "is required to be 0.0, but is "
                << pitchShiftSchedule[0].first;
            mexErrMsgTxt(errMsg.str().c_str());
        }

        // Check that all the time points are monotonically increasing.
        for (size_t i = 1; i < pitchShiftSchedule.size(); ++i) {
            if (pitchShiftSchedule[i].first <= pitchShiftSchedule[i - 1].first) {
                std::ostringstream errMsg;
                errMsg << "The time points in the pitch-shift schedule are "
                    << "required to be monotonically increasing, but the value at "
                    << "index " << i << " (" << pitchShiftSchedule[i].first
                    << ") is <= the value at index " << (i - 1) << " ("
                    << pitchShiftSchedule[i - 1].first << ").";
                mexErrMsgTxt(errMsg.str().c_str());
            }
        }

        // Check that pitchShiftRatio is a positive integer and is within bound.
        for (size_t i = 0; i < pitchShiftSchedule.size(); ++i) {
            const std::pair<dtype, dtype> shiftSpec = pitchShiftSchedule[i];
            if (shiftSpec.second <= 0.0) {
                std::ostringstream errMsg;
                errMsg << "Each pitchShiftRatio is required to be a positive "
                    << "number, but is " << shiftSpec.second << " at index "
                    << i << " of " << pitchShiftSchedule.size();
                mexErrMsgTxt(errMsg.str().c_str());
            }
            if (shiftSpec.second > maxPitchShiftRatio) {
                std::ostringstream errMsg;
                errMsg << "Maximum pitch shift ratio (" << maxPitchShiftRatio
                    << ") is exceeded: " << shiftSpec.second;
                mexErrMsgTxt(errMsg.str().c_str());
            }
        }
    }

    void TimeDomainShifter::checkFrameLen(const int len) const {
        if (len != frameLen) {
            std::ostringstream errMsg;
            errMsg << "len (" << len << ") does not match pre-specified "
                << "frameLen (" << frameLen << ").";
            mexErrMsgTxt(errMsg.str().c_str());
        }
    }


    void TimeDomainShifter::smoothScratch(int begin, int end) {
        for (int i = begin; i < end; ++i) {
            const dtype xMinus2 = scratchBuf[getScratchIndex(i - 2)];
            const dtype xMinus1 = scratchBuf[getScratchIndex(i - 1)];
            const dtype x = scratchBuf[getScratchIndex(i)];
            const dtype xPlus1 = scratchBuf[getScratchIndex(i + 1)];
            const dtype xPlus2 = scratchBuf[getScratchIndex(i + 2)];
            scratchBuf[getScratchIndex(i)] =
                0.1 * xMinus2 + 0.2 * xMinus1 + 0.4 * x + 0.2 * xPlus1
                + 0.1 * xPlus2;
        }
    }

    int TimeDomainShifter::getScratchIndex(int i) {
        while (i < 0) {
            i += scratchLen;
        }
        return i % scratchLen;
    }

}  // namespace audapter