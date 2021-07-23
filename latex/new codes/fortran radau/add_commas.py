
tolerances = [6, 12]

for tolerance in tolerances:
    f = open(f"radau_state_tol_{tolerance}.txt", "r")
    lines = f.read()
    print(lines)
    lines = lines.splitlines()

    res = ""
    for line in lines:
        arr = line.split()
        #print(arr)
        for elem in arr:
            res += elem + ", "
        res += "\n"
    f.close()

    f = open(f"radau_state_tol_{tolerance}.csv", "w")
    f.write(res)
