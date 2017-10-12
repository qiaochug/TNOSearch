# Search for Planet Nine
  Welcome to the Planet Nine research repository! The search uses data from the [Dark Energy Survey](https://www.darkenergysurvey.org/) and focuses on transneptunian objects that are ~200AU from earth. Here is some of the code I wrote for the project during summer 2017.
  
## Projects
- trjectory simulation
- position prediction
- candidate evaluation website
- merge triplets algorithm
- analysis of orbit fitter behavior
- Machine Learning orbit classifier

### trajectory simulation
- __trajectory described in angles__: [position.py](position.py)
- __animation of trajectory__: [Asteroid_animated.py](Asteroid_animated.py)
- __generated demos__: [newphiAtane0.4.mp4](newphiAtane0.4.mp4), [ObesrvingFromEquator0e0.4i.mp4](ObesrvingFromEquator0e0.4i.mp4), [ObserveFromNorthPoleFacingRighthandLookingUp0.4e.mp4](ObserveFromNorthPoleFacingRighthandLookingUp0.4e.mp4)
                   
### position prediction
- __prediction with known perilelion__: [predict_by_perilelion.py](predict_by_perilelion.py)
- __prediction based on earth's orbit motion__: [GammaTPlotwStatTNO.py](GammaTPlotwStatTNO.py), [GammaTPlotwStatTNOExFaster.py](GammaTPlotwStatTNOExFaster.py)
- __prediction optimization using inclination__: [InclinationFromTwoObs.py](InclinationFromTwoObs.py)
- __prediction optimization using speed__: [fakeSpeed](speedChange/fakeSpeed)
- __generated demos__: [GammaAgnTEx1--gammaBeta.mp4](GammaAgnTEx1--gammaBeta.mp4), [GammaAgnTEx2--lambdaBeta.mp4](GammaAgnTEx2--lambdaBeta.mp4)

### candidate evaluation website
- __organization of image files__: [batchUntar](linuxScripts/batchUntar), [csvGeneration](linuxScripts/csvGeneration)
- __web app that displays candidate images and records result in cvs__: [web2.html](Website/web2.html)
- __website demo__: [displayexample.png](Website/displayexample.png)

### merge triplets algorithm
- __orbit merging by growing the tree__: [mergeTrips2.py](mergeTriplets/mergeTrips2.py)
- __criteria for new orbit__: [worthadding2.py](mergeTriplets/worthadding2.py), [worthadding3.py](mergeTriplets/worthadding3.py)

### analysis of orbit fitter behavior
- __covariance checking__: [covarDistribution.py](covarCheck/covarDistribution.py)
- __chi-square checking__: [chiSqDistribution.py](chisquareCheck/chiSqDistribution.py)
- __visualization in search of failure pattern__: [chiSqHist.py](chisquareCheck/chiSqHist.py)

### Machine Learning orbit classifier
- __queries into database and features processing__: [datasetprep](ML/datasetprep)
- __Matlab Logistic Regression Model__: [logisticRegression](ML/logisticRegression)
- __Python Random Forest Model__: [randomForest](ML/randomForest)
