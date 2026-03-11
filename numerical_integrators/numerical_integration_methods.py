def forward_euler(f, t, x, h, amod):
    for i in range(1, len(t)):
        x[:, i] = x[:, i-1] + h * f(t[i-1], x[:, i-1], amod)

    return t, x
