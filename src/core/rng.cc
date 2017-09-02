#include "core/rng.h"


#include <typeinfo>


#include "boost/numeric/conversion/cast.hpp"


#include "common_macros.h"


namespace core
{


uint32_t RNG::xsh_rs::apply(uint64_t _input)
{
	const Bitcount_t output_shift = 
		static_cast<Bitcount_t>(_input >> (kInputSize_ - kOpcodeSize)) & kOpcodeMask;
	const uint64_t xorshifted_input = _input ^ (_input >> kInputShift);
	const uint32_t result =
		static_cast<uint32_t>(xorshifted_input >> (kBottomSpare - kOpcodeMask + output_shift));
	return result;
}


uint32_t RNG::xsh_rr::apply(uint64_t _input)
{
	const Bitcount_t output_rotation =
		static_cast<Bitcount_t>(_input >> (kInputSize_ - kOpcodeSize)) & kOpcodeMask;
	const uint64_t xorshifted_input = _input ^ (_input >> kInputShift);
	const uint32_t truncated_input = static_cast<uint32_t>(xorshifted_input >> kBottomSpare);
	return rotate_right(truncated_input, output_rotation);
}


RNG::rxs_m_xs::BoolUint32Pair_t RNG::rxs_m_xs::inverse_step(uint32_t _output, size_t _i)
{
	const uint32_t state = inverse(_output) * kMultiplier + kIncrement
		+ static_cast<uint32_t>(_i * 2);
	const uint32_t result = apply(state);
	const bool is_zero = (state & 3u);
	return { is_zero, result };
}


uint32_t RNG::rxs_m_xs::apply(uint32_t _input)
{
	const Bitcount_t random_shift =
		static_cast<Bitcount_t>(_input >> (32u - kOpcodeSize)) & kOpcodeMask;
	const uint32_t xorshifted_input = xorshift(_input, kOpcodeSize + random_shift);
	const uint32_t multiplied_rxs_input = xorshifted_input * kMcgMultiplier;
	return xorshift(multiplied_rxs_input, kXsShift);
}


uint32_t RNG::rxs_m_xs::inverse(uint32_t _output)
{
	const uint32_t ixs_output = inv_xorshift(_output, kXsShift);
	const uint32_t im_ixs_output = ixs_output * kInvMcgMultiplier;
	const Bitcount_t random_shift =
		static_cast<Bitcount_t>((im_ixs_output >> (32u - kOpcodeSize)) & kOpcodeMask);
	return inv_xorshift(im_ixs_output, kOpcodeSize + random_shift);
}


RNG::RNG(uint64_t _seed) :
	state_{ (_seed + kIncrement_) * kMultiplier_ + kIncrement_ },
	extension_{}
{
	uint32_t xdiff = generator_value_() - generator_value_();
	for (size_t extension_index = 0; extension_index < kDimensionCount_; ++extension_index)
	{
		extension_[extension_index] = generator_value_() ^ xdiff;
	}
}


uint32_t RNG::operator()()
{
	const uint32_t rhs = extension_value_();
	const uint32_t lhs = generator_value_();
	return lhs ^ rhs;
}


uint32_t RNG::generator_value_()
{
	const uint64_t current_state = state_;
	state_ = state_ * kMultiplier_ + kIncrement_;
	return xsh_rr::apply(current_state);
}


uint32_t RNG::extension_value_()
{
	const size_t extension_index = static_cast<size_t>(state_ & kDimensionMask_);
	const bool does_tick = ((state_ & kTickMask_) == 0u);
	if (does_tick)
		advance_extension_();
	return extension_[extension_index];
}


uint32_t RNG::xorshift(uint32_t _input, Bitcount_t _shift)
{
	return _input ^ (_input >> _shift);
}


uint32_t RNG::inv_xorshift(uint32_t _output, Bitcount_t _shift)
{
	return __inv_xorshift_(_output, 32u, _shift);
}


uint32_t RNG::rotate_right(uint32_t _input, Bitcount_t _rotation)
{
	constexpr Bitcount_t kMask = 32u - 1u;
	return (_input >> _rotation) | (_input << ((-_rotation) & kMask));
}


uint32_t RNG::__inv_xorshift_(uint32_t _output, Bitcount_t _bitcount, Bitcount_t _shift)
{
	if (2 * _shift >= _bitcount)
		return xorshift(_output, _shift);
	const uint32_t lowmask1 = (1u << (_bitcount - _shift * 2)) - 1u;
	const uint32_t highmask1 = ~lowmask1;
	const uint32_t bottom1 = _output & lowmask1;
	const uint32_t top = xorshift(_output, _shift) & highmask1;
	const uint32_t combined1 = top | bottom1;
	const uint32_t lowmask2 = (1u << (_bitcount - _shift)) - 1u;
	const uint32_t bottom2 = combined1 & lowmask2;
	const uint32_t inv_bottom2 = __inv_xorshift_(bottom2, _bitcount - _shift, _shift);
	return top | (inv_bottom2 & lowmask1);
}


void RNG::advance_extension_()
{
	bool carry = false;
	for (ExtensionArray_t::iterator ecit = extension_.begin();
		 ecit != extension_.end(); ++ecit)
	{
		const size_t i = boost::numeric_cast<size_t>(std::distance(extension_.begin(), ecit));
		uint32_t& ext = *ecit;
		if (carry)
		{
			const rxs_m_xs::BoolUint32Pair_t step_result = rxs_m_xs::inverse_step(ext, i + 1u);
			carry = std::get<0>(step_result);
			ext = std::get<1>(step_result);
		}
		const rxs_m_xs::BoolUint32Pair_t step_result = rxs_m_xs::inverse_step(ext, i + 1u);
		ext = std::get<1>(step_result);
		carry = carry || std::get<0>(step_result);
	}
}


} // namespace core
