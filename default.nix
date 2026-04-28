{ stdenv, cmake, pkg-config, onetbb, matplotplusplus, biu }: stdenv.mkDerivation
{
  name = "ufo";
  src = ./.;
  buildInputs = [ onetbb matplotplusplus biu ];
  nativeBuildInputs = [ cmake pkg-config ];
  doCheck = true;
}
