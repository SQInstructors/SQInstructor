from sqitorsion import SQItorsion, SQIscalar, SQIborel

def sqi_test():
    # sqi = SQItorsion()
    # sqi = SQIscalar()
    sqi = SQIborel()
    sqi.keygen()
    print("Keygen done")
    sqi.commitment()
    print("Commitment done")
    # M = sqi.challenge()
    m = sqi.challenge(debug=True)
    print("Challenge done")
    # phi = sqi.response(M)
    phi = sqi.response(m)
    print("Equivalent ideal done")
    # out = sqi.verify(M, phi)
    out = sqi.verify(m, phi)
    print(f"Verification output: {out}")
    print("Verification done")

if __name__=="__main__":
    print("Starting sqi_test")
    sqi_test()
