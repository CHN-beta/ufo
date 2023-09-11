{
  inputs.nixos.url = "github:CHN-beta/nixos";

  outputs = inputs: let pkgs = inputs.nixos.nixosConfigurations.pc.pkgs; in
  {
    devShell.x86_64-linux = pkgs.mkShell.override { stdenv = pkgs.genericPackages.gcc13Stdenv; }
    {
      buildInputs = with pkgs;
        [ yaml-cpp eigen fmt (localPackages.concurrencpp.override { stdenv = genericPackages.gcc13Stdenv; }) ];
      nativeBuildInputs = with pkgs; [ gdb ];
    };
  };
}
