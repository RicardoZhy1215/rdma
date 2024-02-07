#include <stdio.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
/*
*/


int main() {
    pid_t pid = fork();

    if (pid == 0) {
        // Child process
        printf("Child process is running\n");
        printf("Child process is running\n");
        printf("Child process is running\n");

        // Child process logic...
    } else if (pid > 0) {
        // Parent process
         // Parent process waits for child process to finish
        
        printf("Child has finished\n");
        wait(NULL);

    } else {
        // fork failed
        fprintf(stderr, "fork failed\n");
    }

    return 0;
}
