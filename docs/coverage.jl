# Only run coverage from linux release build on travis.
get(ENV, "TRAVIS_OS_NAME", "")       == "linux"   || exit()
get(ENV, "TRAVIS_JULIA_VERSION", "") == "release" || exit()

using Coverage

Codecov.submit(Codecov.process_folder())
