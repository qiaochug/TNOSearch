# Search for Planet Nine

## Projects
- trjectory simulation
- position prediction
- candidate evaluation website
- merge triplets algorithm
- analysis of orbit fitter behavior
- Machine Learning orbit classifier

### trajectory simulation
- __trajectory described in angles__: position.py
- __animation of trajectory__: Asteroid_animated.py
- __generated demos__: newphiAtane0.4.mp4
                   Observe_from_north_pole,facing_righthand,looking_up,x =y,y=x,e=0.4.mp4
                   obesrving_from_equator,x=y,y=z,e=0,i= 0.4.mp4
                   
### position prediction
- __prediction with known perilelion__: predict_by_perilelion.py
- __prediction based on earth's orbit motion__: GammaTPlotwStatTNO.py, GammaTPlotwStatTNOExFaster.py
- __prediction optimization using inclination__: InclinationFromTwoObs.py
- __prediction optimization using speed__: speedChange/fakeSpeed
- __generated demos__: GammaAgnTEx1--gammaBeta.mp4, GammaAgnTEx2--lambdaBeta.mp4

### candidate evaluation website
- __organization of image files__: linuxScripts/batchUntar, linuxScripts/csvGeneration
- __web app that displays candidate images and records result in cvs__: Website/web2.html
- __website demo__: Website/displayexample.png

### merge triplets algorithm
- __orbit merging by growing the tree__: mergeTriplets/mergeTrips2.py
- __criteria for new orbit__: mergeTriplets/worthadding2.py, mergeTriplets/worthadding3.py

### analysis of orbit fitter behavior
- __covariance checking__: covarCheck/covarDistribution.py
- __chi-square checking__: chisquareCheck/chiSqDistribution.py
- __visualization in search of failure pattern__: chisquareCheck/chiSqHist.py

### Machine Learning orbit classifier
- __queries into database and features processing__: ML/datasetprep
- __Matlab Logistic Regression Model__: ML/logisticRegression
- __Python Random Forest Model__: ML/randomForest
