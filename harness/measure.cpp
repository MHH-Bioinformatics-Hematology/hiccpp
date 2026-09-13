// hicfilecpp-measure: runs a command and reports its CPU time and peak RSS.
//
//   hicfilecpp-measure OUT.json TIMEOUT_SECONDS -- COMMAND [ARGS...]
//
// Linux carries a process's RSS high-water mark across exec, so a child forked
// from a large process (the Python runner) would report the runner's RSS. This
// wrapper is small; the command is forked from it, so the reported peak RSS is
// the command's own. CPU time includes the command's waited-for children.

#include <signal.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <unistd.h>

#include <cerrno>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {

pid_t child = -1;

void on_alarm(int) {
    if (child > 0) {
        kill(child, SIGKILL);
    }
}

void forward(int signal_number) {
    if (child > 0) {
        kill(child, signal_number);
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 5 || std::strcmp(argv[3], "--") != 0) {
        std::fprintf(stderr, "usage: hicfilecpp-measure OUT.json TIMEOUT_SECONDS -- COMMAND [ARGS...]\n");
        return 2;
    }
    const char* out_path = argv[1];
    const long timeout = std::strtol(argv[2], nullptr, 10);
    const auto start = std::chrono::steady_clock::now();
    child = fork();
    if (child < 0) {
        std::perror("fork");
        return 2;
    }
    if (child == 0) {
        execvp(argv[4], argv + 4);
        std::perror("execvp");
        _exit(127);
    }
    signal(SIGALRM, on_alarm);
    signal(SIGTERM, forward);
    signal(SIGINT, forward);
    if (timeout > 0) {
        alarm(static_cast<unsigned>(timeout));
    }
    int status = 0;
    rusage usage{};
    while (wait4(child, &status, 0, &usage) < 0) {
        if (errno != EINTR) {
            std::perror("wait4");
            return 2;
        }
    }
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    const int code = WIFEXITED(status) ? WEXITSTATUS(status) : -WTERMSIG(status);
    const double cpu = static_cast<double>(usage.ru_utime.tv_sec + usage.ru_stime.tv_sec) +
                       static_cast<double>(usage.ru_utime.tv_usec + usage.ru_stime.tv_usec) / 1e6;
    std::FILE* out = std::fopen(out_path, "w");
    if (out == nullptr) {
        std::perror("fopen");
        return 2;
    }
    std::fprintf(out, "{\"exit\": %d, \"cpu_s\": %.6f, \"rss_mb\": %.3f, \"wall_s\": %.6f}\n", code, cpu,
                 static_cast<double>(usage.ru_maxrss) / 1024.0, wall);
    std::fclose(out);
    return 0;
}
