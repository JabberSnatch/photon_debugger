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
		uint64_t cdit = 0;
		uint64_t candidate_divisor = 0;
		while ((cdit < i) && (candidate_divisor < (candidate_prime / 2)))
		{
			candidate_divisor = primes[cdit];
			if (candidate_prime % candidate_divisor == 0)
			{
				++candidate_prime;
				cdit = 0;
			}
			else
				++cdit;
		}
		primes[i] = candidate_prime;
		++candidate_prime;
	}
	return primes;
}

template <uint64_t N, uint64_t FirstIndex>
static constexpr std::array<uint64_t, N> get_primes(
	std::array<uint64_t, FirstIndex> first_primes)
{
	using PrimeArray_t = std::array<uint64_t, N>;
	PrimeArray_t primes{};
	uint64_t candidate_prime = first_primes[FirstIndex - 1] + 1u;
	for (uint64_t i = 0; i < N; ++i)
	{
		const uint64_t total_primes_count = i + FirstIndex - 1;
		uint64_t cdit = 0;
		uint64_t candidate_divisor = 0;
		while ((cdit < total_primes_count) && (candidate_divisor < (candidate_prime / 2)))
		{
			candidate_divisor = cdit < FirstIndex ? first_primes[cdit] : primes[cdit - FirstIndex];
			if (candidate_prime % candidate_divisor == 0)
			{
				++candidate_prime;
				cdit = 0;
			}
			else
				++cdit;
		}
		primes[i] = candidate_prime;
		++candidate_prime;
	}
	return primes;
}


template <uint64_t i>
static constexpr uint64_t prime()
{
	uint64_t candidate_prime = prime<i - 1>() + 1;
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
	uint64_t candidate_prime = previous.back() + 1u;
	//PreviousArray::const_iterator pcit = previous.cbegin();
	uint64_t pcit = 0u;
	uint64_t candidate_divisor = 0u;
	while (candidate_divisor <= (candidate_prime / 2))
	{
		candidate_divisor = previous[pcit];
		if (candidate_prime % candidate_divisor == 0u)
		{
			++candidate_prime;
			pcit = 0u;
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

struct prime_rec
{
private:
	template <uint64_t i>
	struct helper
	{
		template <uint64_t N, typename ...values>
		static constexpr auto upward_primes(uint64_t previous, values ...done)
		{
			if constexpr (i < N)
			{
				std::array<uint64_t, sizeof...(done)+1> computed_primes{done..., previous};
				uint64_t candidate_prime = previous + 1u;
				uint64_t candidate_divisor = 0u;
				uint64_t cpit = 0u;
				while((cpit < i) && (candidate_divisor < (candidate_prime / 2)))
				{
					candidate_divisor = computed_primes[cpit];
					if (candidate_prime % candidate_divisor == 0)
					{
						++candidate_prime;
						cpit = 0;
					}
					else
						++cpit;
				}
				return helper<i + 1>::upward_primes<N>(candidate_prime, done..., previous);
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
			return helper<1>::upward_primes<N>(prime<0>());
		}
	};

public:
	template <uint64_t i>
	static constexpr auto get_primes()
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
	while (candidate_divisor <= (candidate_prime / 2))
	{
		candidate_divisor = primes[cdit];
		if (candidate_prime % candidate_divisor == 0)
		{
			++candidate_prime;
			cdit = 0;
		}
		else
			++cdit;
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
