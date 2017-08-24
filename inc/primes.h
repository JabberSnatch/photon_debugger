#pragma once
#ifndef __YS_PRIMES_HPP__
#define __YS_PRIMES_HPP__

#include <array>


namespace algo
{


template <uint64_t N>
static constexpr std::array<uint64_t, N> get_primes()
{
	std::array<uint64_t, N> primes{};
	uint64_t candidate_prime = 2;
	for (uint64_t i = 0; i < N; ++i)
	{
		uint64_t candidate_divisor = 1;
		for (;;)
		{
			++candidate_divisor;
			if (candidate_divisor > candidate_prime / 2)
			{
				primes[i] = candidate_prime;
				++candidate_prime;
				break;
			}
			if (candidate_prime % candidate_divisor == 0)
			{
				++candidate_prime;
				candidate_divisor = 1;
			}
		}
	}
	return primes;
}

struct primes
{
	template <uint64_t i>
	static constexpr uint64_t prime()
	{
		uint64_t candidate_prime = primes::prime<i - 1>() + 1;
		uint64_t candidate_divisor = 1;
		for (;;)
		{
			++candidate_divisor;
			if (candidate_divisor > candidate_prime / 2)
			{
				return candidate_prime;
			}
			if (candidate_prime % candidate_divisor == 0)
			{
				++candidate_prime;
				candidate_divisor = 1;
			}
		}
	}

	template <>
	static constexpr uint64_t prime<0>()
	{
		return 2;
	}

private:
	template <uint64_t i>
	struct helper
	{
		template <uint64_t N, typename ...values>
		static constexpr auto upward_primes(uint64_t previous, values ...done)
		{
			if constexpr (i < N)
			{
				uint64_t candidate_prime = previous + 1;
				uint64_t candidate_divisor = 1;
				for (;;)
				{
					++candidate_divisor;
					if (candidate_divisor > candidate_prime / 2)
					{
						break;
					}
					if (candidate_prime % candidate_divisor == 0)
					{
						++candidate_prime;
						candidate_divisor = 1;
					}
				}
				return primes::helper<i + 1>::upward_primes<N>(candidate_prime, done..., previous);
			}
			else
			{
				return std::array<uint64_t, N>{done..., previous};
			}
		}
	};

	template <>
	struct helper<0>
	{
		template <uint64_t N, typename ...values>
		static constexpr auto upward_primes(uint64_t previous, values ...done)
		{
			return primes::helper<1>::upward_primes<N>(primes::prime<0>());
		}
	};

public:
	template <uint64_t i>
	static constexpr auto get_array()
	{
		return helper<0>::upward_primes<i>(0);
	}
};


} // namespace algo


#endif // __YS_PRIMES_HPP__
