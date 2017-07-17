cd ../src
  make reset
  make coverage
  cd -

bash ./test_all.sh
lcov --capture --directory ../src --output-file ../doc/coverage/cov.inf --gcov-tool gcov-7
genhtml ../doc/coverage/cov.inf --output-directory ../doc/coverage/lcov

cd ../src
  make reset
  cd -

open ../doc/coverage/lcov/index.html
