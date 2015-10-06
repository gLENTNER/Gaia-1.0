// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Library/Interpolation.cc
//
// General interpolation objects. Interp1D and Interp2D provide
// Linear and Bilinear Interpolation respectively.

#include <vector>
#include <iostream>
#include <algorithm>

#include <Interpolate.hpp>

namespace Gaia {

namespace Interpolate {

template<class T>
Linear<T>::Linear(const std::vector<T> &x_, const std::vector<T> &y_){

    if ( x_.empty() || y_.empty() )
    throw InterpError("From Linear::Linear(), one or both of the "
    "input vectors was empty!");

    if ( x_.size() != y_.size() )
    throw InterpError("From Linear::Linear(), the input vectors must "
    "be the same length!");

    for (std::size_t i = 1; i < x_.size(); i++) if (x_[i] < x_[i-1])
    throw InterpError("From Linear::Linear(), the first vector is not "
    "in ascending order!");

    x = x_;
    y = y_;

    m.resize(x.size());

    // solve all changes
    for (std::size_t i = 1; i < x.size(); i++)
        m[i] = (y[i] - y[i-1]) / (x[i] - x[i-1]);
}

template<class T>
std::vector<T> Linear<T>::Interpolate(const std::vector<T> &x_){

    if (x_.empty()) throw InterpError("From Linear::Interpolate(), "
    "your new `x` vector is empty!");

    for (std::size_t i = 1; i < x_.size(); i++) if (x_[i] < x_[i-1])
    throw InterpError("From Linear::Interpolate(), your new `x` vector "
    "is not in ascending order!");

    if ( x_[0] < x[0] || x_[ x_.size() - 1] > x[ x.size() - 1] )
    throw InterpError("From Linear::Interpolate(), your new `x` "
    "vector spreads outside the original domain!");

    std::vector<T> result(x_.size());

    for (std::size_t i = 0; i < x_.size(); i++)
        result[i] = Interpolate(x_[i]);

    return result;
}

template<class T>
T Linear<T>::Interpolate(const T &x_){

    // find appropriate interval
    std::size_t i = std::upper_bound(x.begin(), x.end(), x_) - x.begin();

    return m[i] * (x_ - x[i-1]) + y[i-1];
}

template<class T>
BiLinear<T>::BiLinear( const std::vector<T> &x_, const std::vector<T> &y_,
    const std::vector< std::vector<T> > &z_ ){

    //
    // Retain and check valid status of input vectors
    //

    x = x_;
    y = y_;
    z = z_;

    if ( x_.empty() || y_.empty() || z_.empty() )
    throw InterpError("From BiLinear::BiLinear(), one or more of the "
    "input vectors were empty!");

    for (std::size_t i = 1; i < x.size(); i++) if ( x[i] < x[i-1] )
    throw InterpError("From BiLinear::BiLinear(), the `x` vector was "
    "not in ascending order!");

    for (std::size_t i = 1; i < y.size(); i++) if ( y[i] < y[i-1] )
    throw InterpError("From BiLinear::BiLinear(), the `y` vector was "
    "not in ascending order!");

    for (std::size_t i = 1; i < z.size(); i++)
    if ( z[i].size() != z[i-1].size() )
    throw InterpError("From BiLinear::BiLinear(), not all rows in `z` "
    "had an equal length!");

    if ( y.size() != z.size() )
    throw InterpError("From BiLinear::BiLinear(), size of `y` should "
    "equal the rows in `z`!");

    if ( x.size() != z[0].size() )
    throw InterpError("From BiLinear::BiLinear(), the size of `x` "
    "should equal the columns in `z`!");

    // build Linear interpolation objects for each row
    for (std::size_t i = 0; i < z.size(); i++){

        Linear<T> row(x, z[i]);
        partial.push_back(row);
    }
}

template<class T>
std::vector< std::vector<T> > BiLinear<T>::Interpolate(
    const std::vector<T> &x_, const std::vector<T> &y_){

    if (x_.empty() || y_.empty())
    throw InterpError("From BiLinear::Interpolate(), one or more of "
    "your input vectors are empty!");

    for (std::size_t i = 1; i < x_.size(); i++) if (x_[i] < x_[i-1])
    throw InterpError("From BiLinear::Interpolate(), your `x` vector "
    "is not in ascending order!");

    for (std::size_t i = 1; i < y_.size(); i++) if (y_[i] < y_[i-1])
    throw InterpError("From BiLinear::Interpolate(), your `y` vector "
    "is not in ascending order!");

    if ( x_[0] < x[0] || x_[ x_.size() - 1] > x[ x.size() - 1] )
    throw InterpError("From BiLinear::Interpolate(), the new `x` "
    "vector spreads outside the original domain!");

    if ( y_[0] < y[0] || y_[ y_.size() - 1] > y[ y.size() - 1] )
    throw InterpError("From BiLinear::Interpolate(), the new `y` "
    "vector spreads outside the original domain!");

    // initialize `result` matrix
    std::vector< std::vector<T> > result(y_.size(),
        std::vector<T>(x_.size(), 0.0));

    // interpolate at each grid location
    for (std::size_t i = 0; i < y_.size(); i++)
    for (std::size_t j = 0; j < x_.size(); j++)
        result[i][j] = Interpolate(x_[j], y_[i]);;

    return result;
}

template<class T>
T BiLinear<T>::Interpolate(const T &x_, const T &y_){

    // find proper row
    std::size_t i = std::upper_bound(y.begin(), y.end(), y_) - y.begin();

    // solve in `x` first
    T R1 = partial[i-1].Interpolate(x_);
    T R2 = partial[i  ].Interpolate(x_);

    // then solve in `y`
    return R1 + (y_ - y[i-1]) * (R2 - R1) / (y[i] - y[i-1]);
}

template class Linear<float>;
template class Linear<double>;
template class Linear<long double>;

template class BiLinear<float>;
template class BiLinear<double>;
template class BiLinear<long double>;

} // namespace Interpolate

} // namespace Gaia
