#pragma once

#include <cassert>
#include <cerrno>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <sys/prctl.h>
#include <sys/wait.h>
#include <unistd.h>
#include <utility>

namespace chazelle::test {

template <typename Function> bool assertion_aborts(Function&& fn) {
    const pid_t pid = fork();
    if (pid < 0)
        return false;
    if (pid == 0) {
        if (prctl(PR_SET_DUMPABLE, 0) < 0 || std::freopen("/dev/null", "w", stderr) == nullptr)
            std::_Exit(2);
        fn();
        std::_Exit(0);
    }

    int status = 0;
    pid_t result;
    do {
        result = waitpid(pid, &status, 0);
    } while (result < 0 && errno == EINTR);
    return result == pid && WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT;
}

template <typename Function> void require_assertion_abort(Function&& fn) {
    [[maybe_unused]] const bool aborted = assertion_aborts(std::forward<Function>(fn));
    assert(aborted);
}

}
