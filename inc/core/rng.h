#pragma once
#ifndef __YS_RNG_HPP__
#define __YS_RNG_HPP__

#include <cstdint>
#include <array>
#include <tuple>


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
		static uint32_t apply(uint64_t _input);
		//
		static constexpr Bitcount_t	kOpcodeSize = 3; // depends on input size
		static constexpr Bitcount_t	kOpcodeMask = (1u << kOpcodeSize) - 1u;
		static constexpr Bitcount_t	kBottomSpare = kSpareSize_ - kOpcodeSize;
		static constexpr Bitcount_t	kInputShift = kOpcodeSize + (32 + kOpcodeMask) / 2;
	};
	// 64/32 xorshift, random rotate
	struct xsh_rr
	{
		static uint32_t apply(uint64_t _input);
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
		static BoolUint32Pair_t inverse_step(uint32_t _output, size_t _i);
		static uint32_t apply(uint32_t _input);
		static uint32_t inverse(uint32_t _output);
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
	uint32_t operator()();
	static uint32_t xorshift(uint32_t _input, Bitcount_t _shift);
	static uint32_t inv_xorshift(uint32_t _output, Bitcount_t _shift);
	static uint32_t rotate_right(uint32_t _input, Bitcount_t _rotation);
private:
	static uint32_t __inv_xorshift_(uint32_t _output, Bitcount_t _bitcount, Bitcount_t _shift);
	uint32_t generator_value_();
	uint32_t extension_value_();
	void advance_extension_();
private:
	uint64_t			state_;
	ExtensionArray_t	extension_;
};


} // namespace core


#endif // __YS_RNG_HPP__
