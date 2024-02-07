#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>

int main() {
    pid_t pid = fork();
    if (pid == 0) {
        // child process
        printf("This is the child process. PID = %d\n", getpid());
        
    } else if (pid > 0) {
        // parent process
        printf("This is the parent process. PID = %d, child PID = %d\n", getpid(), pid);
        //function2
    } else {
        // fork failed
        fprintf(stderr, "fork failed\n");
    }
    return 0;
}
