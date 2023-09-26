{
  inputs.nixos.url = "github:CHN-beta/nixos";

  outputs = inputs:
    let
      pkgs = inputs.nixos.nixosConfigurations.pc.pkgs.unstablePackages;
      localPackages = import "${inputs.nixos}/local/pkgs" { inherit pkgs; inherit (inputs.nixpkgs) lib; };
    in
    {
      devShell.x86_64-linux = pkgs.mkShell.override { stdenv = pkgs.gcc13Stdenv; }
      {
        packages = with pkgs; [ pkg-config cmake ninja ];
        buildInputs = with pkgs;
        [
          yaml-cpp eigen fmt localPackages.concurrencpp highfive tbb_2021_8.dev localPackages.matplotplusplus
          localPackages.zpp-bits
        ];
        hardeningDisable = [ "all" ];
        # NIX_DEBUG = "1";
      };
    };
}
