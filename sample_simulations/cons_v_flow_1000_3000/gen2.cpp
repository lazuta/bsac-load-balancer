#include <cstdlib>
#include <cstdio>
#include <ctime>

int main() {
	srand(time(NULL));
	printf("1000 1000\n");
	for(int i = 0; i < 1000; ++i) {
		printf("%d %d\n", 1 + 4999 * (i == 499), 1500 * ((i == 0)) * 1000);
	}
	for(int i = 0; i < 1000; ++i) {
		printf("%d %d %d\n", i + 1, ((i + 1) % 1000) + 1, 100);
	}
}