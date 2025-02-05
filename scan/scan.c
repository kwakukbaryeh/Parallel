void exclusive_scan_iterative(int* data, int* end)
{
    int N = end - data;

    // upsweep phase.
    for (int twod = 1; twod < N; twod*=2)
    {
        int twod1 = twod*2;
        parallel_for (int i = 0; i < N; i += twod1)
            data[i+twod1-1] += data[i+twod-1];
    }

    data[N-1] = 0;

    // downsweep phase.
    for (int twod = N/2; twod >= 1; twod /= 2)
    {
        int twod1 = twod*2;
        parallel_for (int i = 0; i < N; i += twod1)
        {
        int t = data[i+twod-1];
        data[i+twod-1] = data[i+twod1-1];
        // change twod1 below to twod to reverse prefix sum.
        data[i+twod1-1] += t;
        }
    }
}