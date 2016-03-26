#include <cstdio>

int main() {
	freopen("sample.out", "r", stdin);
	freopen("sample2.out", "w", stdout);

	printf("[");
	for(int i = 0; i < 100; ++i) {
		int n;
		scanf("%d", &n);
		printf("%d, ", n);
	}
	printf("]\n");
	printf("[");
	for(int i = 0; i < 100; ++i) {
		int n;
		scanf("%d", &n);
		printf("%d, ", n);
	}
	printf("]\n");
	return 0;
}