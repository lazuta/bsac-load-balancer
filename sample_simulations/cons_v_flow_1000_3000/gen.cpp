#include <cstdlib>
#include <cstdio>
#include <ctime>

int main() {
	srand(time(NULL));
	printf("1000 3000\n");
	for(int i = 0; i < 1000; ++i) {
		printf("%d %d\n", 1 + 4999 * (i == 499), 1500 * ((i == 0)) * 1000);
	}
	for(int i = 0; i < 1000; ++i) {
		printf("%d %d %d\n", i + 1, ((i + 1) % 1000) + 1, 100);
		for(int j = 0; j < 2; ++j) {
        	int l = rand() % 1000 + 1;
			int r = rand() % 1000 + 1;
			while(l == r) r = rand() % 1000 + 1;
			
			printf("%d %d %d\n", r, l, 1 + (rand() % 5));
    	}
    }
}