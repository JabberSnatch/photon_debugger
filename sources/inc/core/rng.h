#pragma once
#ifndef __YS_RNG_HPP__
#define __YS_RNG_HPP__

#include <cstdint>
#include <array>
#include <limits>
#include <tuple>

#include "maths/maths.h"


namespace core
{


// PCG implementation.
// Constants are found on https://github.com/imneme/pcg-cpp
// Reference : http://www.pcg-random.org/
// 2-dimensional 32 bit generator
class RNG
{
private:
	using Bitcount_t = uint8_t;
	//
	static constexpr uint64_t	kMultiplier_ = 6364136223846793005ull;
	static constexpr uint64_t	kIncrement_ = 1442695040888963407ull;
	//
	static constexpr Bitcount_t	kDimensionPow2_ = 1u;
	static constexpr Bitcount_t	kAdvancePow2_ = 16u;
	static constexpr size_t		kDimensionCount_ = (1u << kDimensionPow2_);
	static constexpr uint64_t	kDimensionMask_ = kDimensionCount_ - 1u; // two-dimensional pcg
	//
	static constexpr Bitcount_t	kOutputSize_ = 32u;
	static constexpr Bitcount_t	kInputSize_ = 64u;
	static constexpr Bitcount_t	kSpareSize_ = kInputSize_ - kOutputSize_;
	//
	// may_tick = true
	static constexpr size_t		kTickShift_ = 64u - kAdvancePow2_;
	static constexpr uint64_t	kTickMask_ = (1u << kAdvancePow2_) - 1u;
	// may_tock = false
private:
	using ExtensionArray_t = std::array<uint32_t, kDimensionCount_>;
private:
	// 64/32 xorshift, random shift
	struct xsh_rs
	{
		static uint32_t Apply(uint64_t _input);
		//
		static constexpr Bitcount_t	kOpcodeSize = 3; // depends on input size
		static constexpr Bitcount_t	kOpcodeMask = (1u << kOpcodeSize) - 1u;
		static constexpr Bitcount_t	kBottomSpare = kSpareSize_ - kOpcodeSize;
		static constexpr Bitcount_t	kInputShift = kOpcodeSize + (32 + kOpcodeMask) / 2;
	};
	// 64/32 xorshift, random rotate
	struct xsh_rr
	{
		static uint32_t Apply(uint64_t _input);
		//
		static constexpr Bitcount_t	kOpcodeSize = 5; // depends on input size
		static constexpr Bitcount_t	kOpcodeMask = (1u << kOpcodeSize) - 1u;
		static constexpr Bitcount_t	kBottomSpare = kSpareSize_ - kOpcodeSize;
		static constexpr Bitcount_t	kInputShift = (kOpcodeSize + 32u) / 2;
	};
	//
	// 32/32 random xorshift, mcg multiply, xorshift
	struct rxs_m_xs
	{
		using BoolUint32Pair_t = std::tuple<bool, uint32_t>;
		static BoolUint32Pair_t InverseStep(uint32_t _output, size_t _i);
		static uint32_t Apply(uint32_t _input);
		static uint32_t Inverse(uint32_t _output);
		//
		static constexpr uint32_t kIncrement = 2891336453u;
		static constexpr uint32_t kMultiplier = 747796405u;
		//
		static constexpr uint32_t kMcgMultiplier = 277803737u;
		static constexpr uint32_t kInvMcgMultiplier = 2897767785u;
		//
		static constexpr Bitcount_t kOpcodeSize = 4; // depends on input size
		static constexpr Bitcount_t kOpcodeMask = (1u << kOpcodeSize) - 1u;
		//
		static constexpr Bitcount_t kXsShift = (2u * 32u + 2u) / 3u;
	};
public:
	explicit RNG(uint64_t _seed);
	uint32_t Get32b();
	uint32_t Get32b(uint32_t _max);
	uint64_t Get64b();
	uint64_t Get64b(uint64_t _max);
	float GetSingle();
	double GetDouble();
	maths::Decimal GetDecimal() { return GetFloat<maths::Decimal>(); }
	//
	static uint32_t Xorshift(uint32_t _input, Bitcount_t _shift);
	static uint32_t InvXorshift(uint32_t _output, Bitcount_t _shift);
	static uint32_t RotateRight(uint32_t _input, Bitcount_t _rotation);
private:
	template <typename T> T GetFloat();
	template <> float GetFloat() { return GetSingle(); }
	template <> double GetFloat() { return GetDouble(); }
	static uint32_t __InvXorshift_(uint32_t _output, Bitcount_t _bitcount, Bitcount_t _shift);
	uint32_t GeneratorValue_();
	uint32_t ExtensionValue_();
	void AdvanceExtension_();
private:
	uint64_t			state_;
	ExtensionArray_t	extension_;
};


} // namespace core


#endif // __YS_RNG_HPP__
