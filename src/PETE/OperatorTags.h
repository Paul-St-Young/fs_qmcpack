// -*- C++ -*-
// ACL:license
// ----------------------------------------------------------------------
// This software and ancillary information (herein called "SOFTWARE")
// called PETE (Portable Expression Template Engine) is
// made available under the terms described here.  The SOFTWARE has been
// approved for release with associated LA-CC Number LA-CC-99-5.
//
// Unless otherwise indicated, this SOFTWARE has been authored by an
// employee or employees of the University of California, operator of the
// Los Alamos National Laboratory under Contract No.  W-7405-ENG-36 with
// the U.S. Department of Energy.  The U.S. Government has rights to use,
// reproduce, and distribute this SOFTWARE. The public may copy, distribute,
// prepare derivative works and publicly display this SOFTWARE without
// charge, provided that this Notice and any statement of authorship are
// reproduced on all copies.  Neither the Government nor the University
// makes any warranty, express or implied, or assumes any liability or
// responsibility for the use of this SOFTWARE.
//
// If SOFTWARE is modified to produce derivative works, such modified
// SOFTWARE should be clearly marked, so as not to confuse it with the
// version available from LANL.
//
// For more information about PETE, send e-mail to pete@acl.lanl.gov,
// or visit the PETE web page at http://www.acl.lanl.gov/pete/.
// ----------------------------------------------------------------------
// ACL:license

#include <cstdlib>
#include <cmath>

#ifndef PETE_PETE_OPERATORTAGS_H
#define PETE_PETE_OPERATORTAGS_H

namespace qmcplusplus
{
///////////////////////////////////////////////////////////////////////////////
//
// WARNING: THIS FILE WAS GENERATED AUTOMATICALLY!
// YOU SHOULD MODIFY THE INPUT FILES INSTEAD OF CHANGING THIS FILE DIRECTLY!
//
// THE FOLLOWING INPUT FILES WERE USED TO MAKE THIS FILE:
//
// MakeOperators
// PeteOps.in
//
///////////////////////////////////////////////////////////////////////////////


struct FnArcCos
{
  PETE_EMPTY_CONSTRUCTORS(FnArcCos)
  template<class T>
  inline typename UnaryReturn<T, FnArcCos >::Type_t
  operator()(const T &a) const
  {
    return (acos(a));
  }
};

struct FnArcSin
{
  PETE_EMPTY_CONSTRUCTORS(FnArcSin)
  template<class T>
  inline typename UnaryReturn<T, FnArcSin >::Type_t
  operator()(const T &a) const
  {
    return (asin(a));
  }
};

struct FnArcTan
{
  PETE_EMPTY_CONSTRUCTORS(FnArcTan)
  template<class T>
  inline typename UnaryReturn<T, FnArcTan >::Type_t
  operator()(const T &a) const
  {
    return (atan(a));
  }
};

struct FnCeil
{
  PETE_EMPTY_CONSTRUCTORS(FnCeil)
  template<class T>
  inline typename UnaryReturn<T, FnCeil >::Type_t
  operator()(const T &a) const
  {
    return (ceil(a));
  }
};

struct FnCos
{
  PETE_EMPTY_CONSTRUCTORS(FnCos)
  template<class T>
  inline typename UnaryReturn<T, FnCos >::Type_t
  operator()(const T &a) const
  {
    return (cos(a));
  }
};

struct FnHypCos
{
  PETE_EMPTY_CONSTRUCTORS(FnHypCos)
  template<class T>
  inline typename UnaryReturn<T, FnHypCos >::Type_t
  operator()(const T &a) const
  {
    return (cosh(a));
  }
};

struct FnExp
{
  PETE_EMPTY_CONSTRUCTORS(FnExp)
  template<class T>
  inline typename UnaryReturn<T, FnExp >::Type_t
  operator()(const T &a) const
  {
    return (exp(a));
  }
};

struct FnFabs
{
  PETE_EMPTY_CONSTRUCTORS(FnFabs)
  template<class T>
  inline typename UnaryReturn<T, FnFabs >::Type_t
  operator()(const T &a) const
  {
    return (std::abs(a));
  }
};

struct FnFloor
{
  PETE_EMPTY_CONSTRUCTORS(FnFloor)
  template<class T>
  inline typename UnaryReturn<T, FnFloor >::Type_t
  operator()(const T &a) const
  {
    return (floor(a));
  }
};

struct FnLog
{
  PETE_EMPTY_CONSTRUCTORS(FnLog)
  template<class T>
  inline typename UnaryReturn<T, FnLog >::Type_t
  operator()(const T &a) const
  {
    return (log(a));
  }
};

struct FnLog10
{
  PETE_EMPTY_CONSTRUCTORS(FnLog10)
  template<class T>
  inline typename UnaryReturn<T, FnLog10 >::Type_t
  operator()(const T &a) const
  {
    return (log10(a));
  }
};

struct FnSin
{
  PETE_EMPTY_CONSTRUCTORS(FnSin)
  template<class T>
  inline typename UnaryReturn<T, FnSin >::Type_t
  operator()(const T &a) const
  {
    return (sin(a));
  }
};

struct FnHypSin
{
  PETE_EMPTY_CONSTRUCTORS(FnHypSin)
  template<class T>
  inline typename UnaryReturn<T, FnHypSin >::Type_t
  operator()(const T &a) const
  {
    return (sinh(a));
  }
};

struct FnSqrt
{
  PETE_EMPTY_CONSTRUCTORS(FnSqrt)
  template<class T>
  inline typename UnaryReturn<T, FnSqrt >::Type_t
  operator()(const T &a) const
  {
    return (sqrt(a));
  }
};

struct FnTan
{
  PETE_EMPTY_CONSTRUCTORS(FnTan)
  template<class T>
  inline typename UnaryReturn<T, FnTan >::Type_t
  operator()(const T &a) const
  {
    return (tan(a));
  }
};

struct FnHypTan
{
  PETE_EMPTY_CONSTRUCTORS(FnHypTan)
  template<class T>
  inline typename UnaryReturn<T, FnHypTan >::Type_t
  operator()(const T &a) const
  {
    return (tanh(a));
  }
};

struct OpUnaryMinus
{
  PETE_EMPTY_CONSTRUCTORS(OpUnaryMinus)
  template<class T>
  inline typename UnaryReturn<T, OpUnaryMinus >::Type_t
  operator()(const T &a) const
  {
    return (-a);
  }
};

struct OpUnaryPlus
{
  PETE_EMPTY_CONSTRUCTORS(OpUnaryPlus)
  template<class T>
  inline typename UnaryReturn<T, OpUnaryPlus >::Type_t
  operator()(const T &a) const
  {
    return (+a);
  }
};

struct OpBitwiseNot
{
  PETE_EMPTY_CONSTRUCTORS(OpBitwiseNot)
  template<class T>
  inline typename UnaryReturn<T, OpBitwiseNot >::Type_t
  operator()(const T &a) const
  {
    return (~a);
  }
};

struct OpIdentity
{
  PETE_EMPTY_CONSTRUCTORS(OpIdentity)
  template<class T>
  inline typename UnaryReturn<T, OpIdentity >::Type_t
  operator()(const T &a) const
  {
    return (a);
  }
};

struct OpNot
{
  PETE_EMPTY_CONSTRUCTORS(OpNot)
  template<class T>
  inline typename UnaryReturn<T, OpNot >::Type_t
  operator()(const T &a) const
  {
    return (!a);
  }
};

template<class T >
struct UnaryReturn<T, OpNot >
{
  typedef bool Type_t;
};

template <class T1>
struct OpCast
{
  PETE_EMPTY_CONSTRUCTORS_TEMPLATE(OpCast, T1)
  template<class T2>
  inline UnaryReturn<T2, OpCast<T1> >
  operator()(const T2 &a) const
  {
    return T1(a);
  }
};

template<class T1, class T2>
struct UnaryReturn<T2, OpCast<T1> >
{
  typedef T1 Type_t;
};

struct OpAdd
{
  PETE_EMPTY_CONSTRUCTORS(OpAdd)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpAdd >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a + b);
  }
};

struct OpSubtract
{
  PETE_EMPTY_CONSTRUCTORS(OpSubtract)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpSubtract >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a - b);
  }
};

struct OpMultiply
{
  PETE_EMPTY_CONSTRUCTORS(OpMultiply)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpMultiply >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a * b);
  }
};

struct OpDivide
{
  PETE_EMPTY_CONSTRUCTORS(OpDivide)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpDivide >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a / b);
  }
};

struct OpMod
{
  PETE_EMPTY_CONSTRUCTORS(OpMod)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpMod >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a % b);
  }
};

struct OpBitwiseAnd
{
  PETE_EMPTY_CONSTRUCTORS(OpBitwiseAnd)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpBitwiseAnd >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a & b);
  }
};

struct OpBitwiseOr
{
  PETE_EMPTY_CONSTRUCTORS(OpBitwiseOr)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpBitwiseOr >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a | b);
  }
};

struct OpBitwiseXor
{
  PETE_EMPTY_CONSTRUCTORS(OpBitwiseXor)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpBitwiseXor >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a ^ b);
  }
};

struct FnLdexp
{
  PETE_EMPTY_CONSTRUCTORS(FnLdexp)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, FnLdexp >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (ldexp(a,b));
  }
};

struct FnPow
{
  PETE_EMPTY_CONSTRUCTORS(FnPow)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, FnPow >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (pow(a,b));
  }
};

struct FnFmod
{
  PETE_EMPTY_CONSTRUCTORS(FnFmod)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, FnFmod >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (fmod(a,b));
  }
};

struct FnArcTan2
{
  PETE_EMPTY_CONSTRUCTORS(FnArcTan2)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, FnArcTan2 >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (atan2(a,b));
  }
};

struct OpLT
{
  PETE_EMPTY_CONSTRUCTORS(OpLT)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpLT >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a < b);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpLT >
{
  typedef bool Type_t;
};

struct OpLE
{
  PETE_EMPTY_CONSTRUCTORS(OpLE)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpLE >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a <= b);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpLE >
{
  typedef bool Type_t;
};

struct OpGT
{
  PETE_EMPTY_CONSTRUCTORS(OpGT)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpGT >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a > b);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpGT >
{
  typedef bool Type_t;
};

struct OpGE
{
  PETE_EMPTY_CONSTRUCTORS(OpGE)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpGE >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a >= b);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpGE >
{
  typedef bool Type_t;
};

struct OpEQ
{
  PETE_EMPTY_CONSTRUCTORS(OpEQ)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpEQ >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a == b);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpEQ >
{
  typedef bool Type_t;
};

struct OpNE
{
  PETE_EMPTY_CONSTRUCTORS(OpNE)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpNE >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a != b);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpNE >
{
  typedef bool Type_t;
};

struct OpAnd
{
  PETE_EMPTY_CONSTRUCTORS(OpAnd)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpAnd >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a && b);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpAnd >
{
  typedef bool Type_t;
};

struct OpOr
{
  PETE_EMPTY_CONSTRUCTORS(OpOr)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpOr >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a || b);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpOr >
{
  typedef bool Type_t;
};

struct OpLeftShift
{
  PETE_EMPTY_CONSTRUCTORS(OpLeftShift)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpLeftShift >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a << b);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpLeftShift >
{
  typedef T1 Type_t;
};

struct OpRightShift
{
  PETE_EMPTY_CONSTRUCTORS(OpRightShift)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpRightShift >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (a >> b);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpRightShift >
{
  typedef T1 Type_t;
};

struct OpAddAssign
{
  PETE_EMPTY_CONSTRUCTORS(OpAddAssign)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpAddAssign >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    (const_cast<T1 &>(a) += b);
    return const_cast<T1 &>(a);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpAddAssign >
{
  typedef T1 &Type_t;
};

struct OpSubtractAssign
{
  PETE_EMPTY_CONSTRUCTORS(OpSubtractAssign)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpSubtractAssign >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    (const_cast<T1 &>(a) -= b);
    return const_cast<T1 &>(a);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpSubtractAssign >
{
  typedef T1 &Type_t;
};

struct OpMultiplyAssign
{
  PETE_EMPTY_CONSTRUCTORS(OpMultiplyAssign)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpMultiplyAssign >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    (const_cast<T1 &>(a) *= b);
    return const_cast<T1 &>(a);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpMultiplyAssign >
{
  typedef T1 &Type_t;
};

struct OpDivideAssign
{
  PETE_EMPTY_CONSTRUCTORS(OpDivideAssign)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpDivideAssign >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    (const_cast<T1 &>(a) /= b);
    return const_cast<T1 &>(a);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpDivideAssign >
{
  typedef T1 &Type_t;
};

struct OpModAssign
{
  PETE_EMPTY_CONSTRUCTORS(OpModAssign)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpModAssign >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    (const_cast<T1 &>(a) %= b);
    return const_cast<T1 &>(a);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpModAssign >
{
  typedef T1 &Type_t;
};

struct OpBitwiseOrAssign
{
  PETE_EMPTY_CONSTRUCTORS(OpBitwiseOrAssign)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpBitwiseOrAssign >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    (const_cast<T1 &>(a) |= b);
    return const_cast<T1 &>(a);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpBitwiseOrAssign >
{
  typedef T1 &Type_t;
};

struct OpBitwiseAndAssign
{
  PETE_EMPTY_CONSTRUCTORS(OpBitwiseAndAssign)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpBitwiseAndAssign >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    (const_cast<T1 &>(a) &= b);
    return const_cast<T1 &>(a);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpBitwiseAndAssign >
{
  typedef T1 &Type_t;
};

struct OpBitwiseXorAssign
{
  PETE_EMPTY_CONSTRUCTORS(OpBitwiseXorAssign)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpBitwiseXorAssign >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    (const_cast<T1 &>(a) ^= b);
    return const_cast<T1 &>(a);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpBitwiseXorAssign >
{
  typedef T1 &Type_t;
};

struct OpLeftShiftAssign
{
  PETE_EMPTY_CONSTRUCTORS(OpLeftShiftAssign)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpLeftShiftAssign >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    (const_cast<T1 &>(a) <<= b);
    return const_cast<T1 &>(a);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpLeftShiftAssign >
{
  typedef T1 &Type_t;
};

struct OpRightShiftAssign
{
  PETE_EMPTY_CONSTRUCTORS(OpRightShiftAssign)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpRightShiftAssign >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    (const_cast<T1 &>(a) >>= b);
    return const_cast<T1 &>(a);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpRightShiftAssign >
{
  typedef T1 &Type_t;
};

struct OpAssign
{
  PETE_EMPTY_CONSTRUCTORS(OpAssign)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpAssign >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    return (const_cast<T1 &>(a) = b);
  }
};

template<class T1, class T2 >
struct BinaryReturn<T1, T2, OpAssign >
{
  typedef T1 &Type_t;
};

struct FnWhere
{
  PETE_EMPTY_CONSTRUCTORS(FnWhere)
  template<class T1, class T2, class T3>
  inline typename TrinaryReturn<T1, T2, T3, FnWhere >
  ::Type_t
  operator()(T1 &a, const T2 &b, const T3 &c) const
  {
    if (a)
      return b;
    else
      return c;
  }
};

}
#endif // PETE_PETE_OPERATORTAGS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile$   $Author: berrill $
// $Revision: 6941 $   $Date: 2016-05-27 09:00:55 -0600 (Fri, 27 May 2016) $
// ----------------------------------------------------------------------
// ACL:rcsinfo
