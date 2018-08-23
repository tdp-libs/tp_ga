#ifndef tp_ga_ArrayRefine_h
#define tp_ga_ArrayRefine_h

#include "tp_ga/Globals.h"

#include "tp_math_utils/Globals.h"

namespace tp_ga
{
namespace detail
{
//##################################################################################################
template<class T>
constexpr auto getSize(T const& arr) noexcept -> decltype (arr.length())
{
  return arr.length();
}

//##################################################################################################
template<class T>
constexpr auto getSize(T const& arr) noexcept -> decltype (arr.size())
{
  return arr.size();
}

//##################################################################################################
template<typename DistanceType,
         typename ParameterContainer,
         typename ScoreFunction>
void refineArray(ParameterContainer& params,
                 ScoreFunction score,
                 size_t iterations)
{
  const size_t poolSize{20};

  struct Pair_lt
  {
    DistanceType score{0};
    ParameterContainer params;
    Pair_lt(DistanceType score_, ParameterContainer params_):
      score(score_),
      params(std::move(params_))
    {

    }
  };

  //Assuming that the params are a std::array or glm::vec this should result in the pool being
  //allocated in contiguous memory to save thrashing the cache in the loop below.
  std::array<Pair_lt, poolSize> scratch{tpMakeArray<Pair_lt, poolSize>(Pair_lt(score(params), params))};
  std::vector<Pair_lt*> pool;
  pool.reserve(scratch.size());
  for(Pair_lt& p : scratch)
    pool.push_back(&p);

  tp_math_utils::RNG rng;

  auto blend = [](tp_math_utils::RNG& rng, auto a, auto b)
  {
    decltype (a) f = rng.randF();
    a = (a*f) + (b*(decltype (a)(1)-f));
    if(rng.randF()<0.4f)
      a += decltype (a)(rng.randF2()*0.00001f);
    if(rng.randF()<0.2f)
      a += decltype (a)(rng.randF2()*0.002f);
    if(rng.randF()<0.01f)
      a += decltype (a)(rng.randF2()*0.2f);
    return a;
  };

  for(size_t i=iterations; i; i--)
  {
    //-- Mutate ------------------------------------------------------------------------------------
    Pair_lt* p = pool.back();
    {
      pool.pop_back();

      std::vector<size_t> indexes = rng.randomIndexes(pool.size(), 2);
      if(indexes.size() == 2)
      {
        Pair_lt* a = pool.at(indexes.at(0));
        Pair_lt* b = pool.at(indexes.at(1));
        size_t cMax=getSize(params);
        for(size_t c=0; c<cMax; c++)
          p->params[c] = blend(rng, a->params[c], b->params[c]);
      }
    }

    //-- Score--------------------------------------------------------------------------------------
    p->score = score(p->params);

    //-- Update pool -------------------------------------------------------------------------------
    {
      bool inserted=false;
      for(size_t j=0; j<pool.size(); j++)
      {
        if(p->score<=pool.at(j)->score)
        {
          tpInsert(pool, j, p);
          inserted=true;
          break;
        }
      }

      if(!inserted)
        pool.push_back(p);
    }
  }

  params = pool.front()->params;
}
}

//##################################################################################################
template<class T>
constexpr size_t getSize(T const& arr) noexcept
{
  return size_t(detail::getSize(arr));
}

//##################################################################################################
template<typename DistanceType,
         typename ParameterContainer,
         typename ObservationContainer,
         typename DistanceFunction>
void refineArray(ParameterContainer& params,
                 const ObservationContainer& obs,
                 DistanceFunction dist,
                 size_t iterations)
{
  auto score = [&](const ParameterContainer& p)
  {
    DistanceType score=DistanceType();
    for(const auto& o : obs)
      score += dist(p, o);
    return score;
  };

  detail::refineArray<DistanceType, ParameterContainer>(params, score, iterations);
}

//##################################################################################################
template<typename DistanceType,
         typename ParameterContainer,
         typename ObservationContainer,
         typename PrepareFunction,
         typename DistanceFunction>
void refineArray(ParameterContainer& params,
                 const ObservationContainer& obs,
                 PrepareFunction prepare,
                 DistanceFunction dist,
                 size_t iterations)
{
  auto score = [&](const ParameterContainer& p)
  {
    DistanceType score=DistanceType();
    const auto pp = prepare(p);
    for(const auto& o : obs)
      score += dist(pp, o);
    return score;
  };

  detail::refineArray<DistanceType, ParameterContainer>(params, score, iterations);
}

}

#endif
