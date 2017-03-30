#ifndef __YS_MATHS_HPP__
#define __YS_MATHS_HPP__


namespace maths
{


constexpr float almost_onef() { return 1.f - std::numeric_limits<float>::epsilon(); }

constexpr float minf(float _lhs, float _rhs) { return (_lhs < _rhs) ? _lhs : _rhs; }
constexpr float maxf(float _lhs, float _rhs) { return (_lhs > _rhs) ? _lhs : _rhs; }
constexpr float clampf(float _v, float _min, float _max) { return minf(maxf(_v, _min), _max); }


} // namespace maths


#endif // __YS_MATHS_HPP__
