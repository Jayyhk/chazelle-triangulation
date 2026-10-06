#include "render.h"
#include "animation/config.h"

#include <cerrno>
#include <cstring>
#include <spawn.h>
#include <stdexcept>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

extern char** environ;

namespace chazelle {

void render_animation(const std::filesystem::path& trace, const std::filesystem::path& video) {
    std::vector<std::string> arguments{CHAZELLE_ANIMATION_PYTHON, CHAZELLE_ANIMATION_SCRIPT,
                                       trace.string(), video.string()};
    std::vector<char*> argv;
    argv.reserve(arguments.size() + 1);
    for (std::string& argument : arguments)
        argv.push_back(argument.data());
    argv.push_back(nullptr);

    posix_spawn_file_actions_t actions;
    int error = posix_spawn_file_actions_init(&actions);
    if (error != 0)
        throw std::runtime_error(std::string("Failed to prepare Manim: ") + std::strerror(error));
    error = posix_spawn_file_actions_adddup2(&actions, STDERR_FILENO, STDOUT_FILENO);
    pid_t process = 0;
    if (error == 0)
        error = posix_spawn(&process, argv.front(), &actions, nullptr, argv.data(), environ);
    posix_spawn_file_actions_destroy(&actions);
    if (error != 0)
        throw std::runtime_error(std::string("Failed to start Manim: ") + std::strerror(error) +
                                 ". Install src/animation/requirements.txt in .venv or configure "
                                 "CHAZELLE_ANIMATION_PYTHON.");
    int status = 0;
    while (waitpid(process, &status, 0) == -1) {
        if (errno != EINTR)
            throw std::runtime_error(std::string("Failed to wait for Manim: ") +
                                     std::strerror(errno));
    }
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
        throw std::runtime_error("Manim rendering failed; the exact trace is saved at " +
                                 trace.string());
}

}
