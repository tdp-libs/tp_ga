#ifndef tp_ga_RANSAC_h
#define tp_ga_RANSAC_h

#include "tp_ga/Globals.h"

#include "tp_math_utils/Globals.h"
#include "tp_utils/DebugUtils.h"
#include "tp_utils/TimeUtils.h"

namespace tp_ga
{
template<typename DistanceType,
         typename ModelType,
         typename ObservationContainer,
         typename FitFunction,
         typename DistanceFunction>
void ransac(ModelType& model,
            const ObservationContainer& obs,
            size_t obsPerSample,
            size_t nSamples,
            FitFunction fit,
            DistanceFunction dist)
{
  tp_math_utils::RNG rng;

  auto bestDistance = dist(model, obs);

  std::vector<size_t> indexes;

  for(auto i=nSamples; i; i--)
  {
    ObservationContainer sample;
    {
      rng.randomIndexes(obs.size(), obsPerSample, indexes);
      sample.reserve(indexes.size());
      for(auto index : indexes)
        sample.push_back(obs.at(index));
    }

    {
      ModelType newModel = model;
      fit(newModel, sample);

      auto distance = dist(newModel, obs);
      if(distance<bestDistance)
      {
        bestDistance = distance;
        model = newModel;
      }
    }
  }
}

}

#endif
