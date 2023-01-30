### setup instructions

```
cmsrel CMSSW_10_6_20
cd CMSSW_10_6_20/src
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
cd PhysicsTools/NanoAODTools/python
git clone https://github.com/bchiarito/cmssw_backend_attoplotting.git fmk_plotting
cd ../../..
scram b
```
