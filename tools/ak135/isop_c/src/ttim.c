#include <stdio.h>
#include <math.h>

#define max 500

main(argc, argv)
	int argc;
	char **argv;
{
	int i, j, k, n, mn[max];
	char *phcd[max];
	float tt[max], dtdd[max], dtdh[max], dddp[max], ts[max], p[max];
	float zs, delta;

	for(i = 0; i < max; i++) phcd[i] = (char *)malloc(10);

	if(tabin("iasp91"))
	{
		printf("Cannot open iasp91.hed and iasp91.tbl\n");
		exit(1);
	}
for(;;)
{
	brnset();

	printf("Source depth (km): ");
	scanf("%e", &zs);
	depset(zs);
	printf("Enter delta: ");
	scanf("%e", &delta);
	trtm(delta, &n, tt, p, dtdd, dtdh, dddp, phcd);

	printf("delta = %6.2f  depth = %6.2f  n = %d\n", delta, zs, n);
	for(i = 0; i < n; i++)
	{
		mn[i] = tt[i]/60.;
		ts[i] = .01*(int)(100.*(tt[i] - mn[i]*60.) + .5);
		if(ts[i] >= 60.)
		{
			mn[i] = mn[i] + 1;
			ts[i] = ts[i] - 60.;
		}
	}
	for(i = 0; i < n; i++)
	{
		printf("%4d %s %16.6e%9.2f%5d%7.2f%11.2e%11.2e%11.2e\n",
			i, phcd[i], p[i], tt[i], mn[i], ts[i],
			dtdd[i], dtdh[i], dddp[i]);
	}
}
}
