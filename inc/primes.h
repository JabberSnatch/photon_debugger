#pragma once
#ifndef __YS_PRIMES_HPP__
#define __YS_PRIMES_HPP__

#include <array>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/transform_view.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/back.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/comparison.hpp>
#include <boost/mpl/arithmetic.hpp>


namespace algo
{


template <uint64_t N>
static constexpr std::array<uint64_t, N> get_primes()
{
	using PrimeArray_t = std::array<uint64_t, N>;
	PrimeArray_t primes{};
	uint64_t candidate_prime = 2;
	for (uint64_t i = 0; i < N; ++i)
	{
		//PrimeArray_t::const_iterator cdit = primes.cbegin();
		uint64_t cdit = 0;
		for (;;)
		{
			//uint64_t candidate_divisor = *cdit;
			uint64_t candidate_divisor = primes[cdit];
			if (candidate_divisor > candidate_prime / 2 || cdit >= i)
			{
				primes[i] = candidate_prime;
				++candidate_prime;
				break;
			}
			if (candidate_prime % candidate_divisor == 0)
			{
				++candidate_prime;
				cdit = 0;
			}
			else
				++cdit;

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

	template <typename ...Values>
	static constexpr uint64_t compute_next_prime(Values ..._previous)
	{
		using PreviousArray = std::array<uint64_t, sizeof...(_previous)>;
		PreviousArray previous{_previous...};
		uint64_t candidate_prime = previous.back() + 1;
		//PreviousArray::const_iterator pcit = previous.cbegin();
		uint64_t pcit = 0;
		for (;;)
		{
			uint64_t candidate_divisor = previous[pcit];
			if (candidate_divisor >= candidate_prime / 2)
				break;
			if (candidate_prime % candidate_divisor == 0)
			{
				++candidate_prime;
				pcit = 0;
			}
			else
				++pcit;
		}
		return candidate_prime;
	}
	template <>
	static constexpr uint64_t compute_next_prime()
	{
		return 2;
	}

	template <uint64_t P>
	struct Prime
	{
		using tag = boost::mpl::integral_c_tag;
		using value_type = uint64_t;
		static constexpr value_type value{ P };
		using type = Prime;
		using next = Prime<value + 1>;
		using prior = Prime<value - 1>;
		constexpr operator uint64_t() const { return (this->value); }
	};


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


#if 0
template <typename Values>
struct ComputeNextPrime
{
	template <typename CV, typename PPI>
	struct Step
	{
		using CandidateValue = CV;
		using PriorPrimeIter = PPI;
		using PriorPrimeValue = boost::mpl::deref<PriorPrimeIter>::type;
		using HalfCandidateValue = boost::mpl::divides<
			CandidateValue, boost::mpl::integral_c<uint64_t, 2>>::type;
		using FoundTest = boost::mpl::greater_equal<PriorPrimeValue, HalfCandidateValue>::type;
		using CandidateModPriorValue = boost::mpl::modulus<CandidateValue, PriorPrimeValue>::type;
		using NextCandidateTest = boost::mpl::equal_to<
			CandidateModPriorValue, boost::mpl::integral_c<uint64_t, 0>>::type;

		using StepDefault = Step<CandidateValue, boost::mpl::next<PriorPrimeIter>::type>;
		using StepNextCandidate = Step<
			boost::mpl::next<CandidateValue>::type, boost::mpl::begin<Values>::type>;
		using StepPicked = boost::mpl::if_<NextCandidateTest, StepNextCandidate, StepDefault>::type;

		using Next = boost::mpl::if_<FoundTest, CandidateValue, StepPicked>::type;
	};

	using BeginCandidateValue = boost::mpl::next<boost::mpl::back<Values>::type>::type;
	using BeginPriorPrimeIter = boost::mpl::begin<Values>::type;
	using Begin = Step<BeginCandidateValue, BeginPriorPrimeIter>;
	/*
	using PreviousArray = std::array<uint64_t, sizeof...(_previous)>;
	PreviousArray previous{ _previous... };
	uint64_t candidate_prime = previous.back() + 1;
	PreviousArray::const_iterator pcit = previous.cbegin();
	for (;;)
	{
		if (*pcit >= candidate_prime / 2)
			break;
		if (candidate_prime % *pcit == 0)
		{
			++candidate_prime;
			pcit = previous.cbegin();
		}
		else
			++pcit;
	}
	return candidate_prime;
	*/
};

struct test
{
	template <uint64_t N>
	using prime = boost::mpl::integral_c<uint64_t, N>;
	static const uint64_t muc = prime<8>::value;
	using prime_vector = boost::mpl::vector_c<uint64_t, 2, 3, 5>;
	using test_compute = ComputeNextPrime<prime_vector>;
	
	struct begin
	{
		using step = test_compute::Begin;
		using candidate_prime = step::CandidateValue;
		using pcit = step::PriorPrimeValue;
		using half_candidate = step::HalfCandidateValue;
		using found_test_result = step::FoundTest;
		using modulus_result = step::CandidateModPriorValue;
	};

	struct next
	{
		using step = begin::step::StepPicked;
		using candidate_prime = step::CandidateValue;
		using pcit = step::PriorPrimeValue;
		using half_candidate = step::HalfCandidateValue;
		using found_test_result = step::FoundTest;
		using modulus_result = step::CandidateModPriorValue;
	};
};
#endif

} // namespace algo


#endif // __YS_PRIMES_HPP__
