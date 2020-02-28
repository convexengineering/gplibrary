export PIP=`which pip`
echo $PIP
python -c "from gpkit.tests.test_repo import test_repos; test_repos(xmloutput=True, ingpkitmodels=True)"
