#ifndef MT7_GAS_PAR_H
#define MT7_GAS_PAR_H

#include <RcppParallel.h>
#include <fenv.h>
#include <cstddef>

// The thread worker shared by the scalar and the joint score-driven filters,
// which run their groups over threads on the fast route.

namespace mt7 {

// below this many groups the opening of the parallel region is not worth
// its cost; measured on the panel benchmarks and held to the plan's
// asymmetry argument (a factor of two either way costs almost nothing)
constexpr int kMinGroupsPar = 8;

// The same shape as distributions7's d7::par_for(), and for the same two
// reasons. The loop is NOINLINE and the sequential branch runs THROUGH the
// worker, so both branches execute one compiled copy rather than two the
// compiler is free to optimize apart. And the worker installs the calling
// thread's floating-point environment before its chunk: R's psigamma,
// bessel_k, pgamma and pbeta return per-thread last bits without it
// (measured 2026-08-21), and while the covered fast route reads only
// digamma and trigamma, which are stable, a filter's recursion accumulates
// over time and would carry any such difference forward rather than
// confining it to one element.
#if defined(__GNUC__) || defined(__clang__)
#define MT7_NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define MT7_NOINLINE __declspec(noinline)
#else
#define MT7_NOINLINE
#endif

template <typename Body>
struct GroupWorker : public RcppParallel::Worker {
    const Body& body;
    fenv_t env;
    explicit GroupWorker(const Body& b) : body(b) { fegetenv(&env); }
    MT7_NOINLINE void operator()(std::size_t begin, std::size_t end) {
        fesetenv(&env);
        for (std::size_t g = begin; g < end; ++g) body(g);
    }
};

} // namespace mt7

#endif
