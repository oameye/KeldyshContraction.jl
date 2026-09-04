using Preferences: set_preferences!
set_preferences!("KeldyshContraction", "dispatch_doctor_mode" => "error")

using KeldyshContraction
using ParallelTestRunner: ParallelTestRunner

testsuite = ParallelTestRunner.find_tests(@__DIR__)
args = ParallelTestRunner.parse_args(ARGS)

if ParallelTestRunner.filter_tests!(testsuite, args)
    delete!(testsuite, "quality/JET")
end

ParallelTestRunner.runtests(KeldyshContraction, args; testsuite)
