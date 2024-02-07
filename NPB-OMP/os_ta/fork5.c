#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>

// int main() {
//     int n = 5; // 

//     for (int i = 0; i < n; i++) {
//         pid_t pid = fork();

//         if (pid < 0) {
//             fprintf(stderr, "fork failed\n");
//             return 1;
//         } else if (pid == 0) {
//             printf("This is child process with PID %d and parent PID %d\n", getpid(), getppid());
//             _exit(0); // 
//         }
      
//     }
//     return 0;
// }

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main() {
    int n = 5; 
    pid_t pid;

    for (int i = 0; i < n; i++) {
        pid = fork();
        if (pid < 0) {
            fprintf(stderr, "fork failed\n");
            exit(1);
        } else if (pid == 0) {
           
            printf("Child process %d: My PID is %d and my parent's PID is %d\n", i+1, getpid(), getppid());
            switch(i) {
                case 0:
                    printf("Child 1 is doing its task\n");
                    break;
                case 1:
                    printf("Child 2 is doing its task\n");
                    break;
                case 2:
                    printf("Child 3 is doing its task\n");
                    break;
                case 3:
                    printf("Child 4 is doing its task\n");
                    break;
                case 4:
                    printf("Child 5 is doing its task\n");
                    break;
            }
            exit(0); 
        }


}

return 0;
}
