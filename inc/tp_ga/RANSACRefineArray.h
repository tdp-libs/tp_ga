#ifndef tp_ga_RANSACRefineArray_h
#define tp_ga_RANSACRefineArray_h

#include "tp_ga/RANSAC.h"
#include "tp_ga/RefineArray.h"

namespace tp_ga
{

//##################################################################################################
template<typename DistanceType,
         typename ParameterContainer,
         typename ObservationContainer,
         typename DistanceFunction>
void ransacRefineArray(ParameterContainer& params,
                       const ObservationContainer& obs,
                       DistanceFunction dist,
                       size_t iterations,
                       size_t obsPerSample,
                       size_t nSamples)
{
  auto fit = [&dist, iterations](ParameterContainer& params, const ObservationContainer& sample)
  {
    refineArray<DistanceType>(params, sample, dist, iterations);
  };

  auto containerDist = [&dist](ParameterContainer& p, const ObservationContainer& obs)
  {
    DistanceType totalDist=DistanceType();
    for(const auto& o : obs)
      totalDist += dist(p, o);
    return totalDist;
  };

  ransac<DistanceType>(params, obs, obsPerSample, nSamples, fit, containerDist);
}

//##################################################################################################
template<typename DistanceType,
         typename ParameterContainer,
         typename ObservationContainer,
         typename PrepareFunction,
         typename DistanceFunction>
void ransacRefineArray(ParameterContainer& params,
                       const ObservationContainer& obs,
                       PrepareFunction prepare,
                       DistanceFunction dist,
                       size_t iterations,
                       size_t obsPerSample,
                       size_t nSamples)
{
  auto fit = [&prepare, &dist, iterations](ParameterContainer& params, const ObservationContainer& sample)
  {
    refineArray<DistanceType>(params, sample, prepare, dist, iterations);
  };

  auto containerDist = [&dist, &prepare](ParameterContainer& p, const ObservationContainer& obs)
  {
    const auto pp = prepare(p);
    DistanceType totalDist=DistanceType();
    for(const auto& o : obs)
      totalDist += dist(pp, o);
    return totalDist;
  };

  ransac<DistanceType>(params, obs, obsPerSample, nSamples, fit, containerDist);
}

//##################################################################################################
template<typename DistanceType,
         typename ParameterContainer,
         typename ObservationContainer,
         typename PrepareFunction,
         typename DistanceFunction>
void ransacRefineArray(ParameterContainer& params,
                       const ObservationContainer& obs,
                       PrepareFunction prepare,
                       DistanceFunction dist,
                       size_t iterations,
                       size_t obsPerSample,
                       size_t nSamples,
                       std::vector<std::array<DistanceType, 2>>& distances)
{
  auto fit = [&prepare, &dist, &distances, iterations](ParameterContainer& params, const ObservationContainer& sample)
  {
    refineArray<DistanceType>(params, sample, prepare, dist, iterations, distances);
  };

  auto containerDist = [&dist, &prepare](ParameterContainer& p, const ObservationContainer& obs)
  {
    const auto pp = prepare(p);
    DistanceType totalDist=DistanceType();
    for(const auto& o : obs)
      totalDist += dist(pp, o);
    return totalDist;
  };

  ransac<DistanceType>(params, obs, obsPerSample, nSamples, fit, containerDist);
}

}

#endif
