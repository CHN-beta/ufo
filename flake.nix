{
  inputs.nixos.url = "github:CHN-beta/nixos";
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs = inputs:
    let
      pkgs = inputs.nixpkgs.legacyPackages.x86_64-linux;
      localPackages = inputs.nixos.nixosConfigurations.pc.pkgs.localPackages;
    in
    {
      devShell.x86_64-linux = pkgs.mkShell.override { stdenv = pkgs.stdenvNoCC; }
      {
        packages = with pkgs; [ xmake gcc13 pkg-config cmake ];
        inputsFrom = with pkgs;
        [
          yaml-cpp eigen fmt localPackages.concurrencpp highfive
          hdf5.dev tbb_2021_8.dev localPackages.matplotplusplus
          localPackages.zpp-bits
        ];
        PKG_CONFIG_PATH = "${pkgs.tbb_2021_8.dev}/lib/pkgconfig:${pkgs.yaml-cpp}/share/pkgconfig";
        yaml-cpp_DIR = "${pkgs.yaml-cpp}/share/cmake/yaml-cpp";
        hardeningDisable = [ "all" ];
        NIX_DEBUG = "1";
      };
    };
}
