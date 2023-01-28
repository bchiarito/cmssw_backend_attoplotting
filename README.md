### setup instructions

```
git clone https://github.com/bchiarito/cmssw_temp_attoframework_backend.git backend_plotting
cd backend_plotting
cmsrel CMSSW_10_6_20
cd CMSSW_10_6_20/src
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
mv ../../hthatAnalysis.py PhysicsTools/NanoAODTools/python/postprocessing/modules/
mv ../../myAnalysis.py PhysicsTools/NanoAODTools/python/postprocessing/modules/
mv ../../myCutflow.py PhysicsTools/NanoAODTools/python/postprocessing/modules/
mv ../../postprocessor.py PhysicsTools/NanoAODTools/python/postprocessing/framework/
scram b
cd ../..
```
