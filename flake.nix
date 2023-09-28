{
  inputs.nixos.url = "github:CHN-beta/nixos";
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs = inputs:
    let
      # pkgs = inputs.nixpkgs.legacyPackages.x86_64-linux;
      pkgs = import inputs.nixpkgs
      {
        localSystem = { system = "x86_64-linux"; gcc = { arch = "alderlake"; tune = "alderlake"; }; };
        config.allowUnfree = true;
      };
      localPackages = import "${inputs.nixos}/local/pkgs" { inherit pkgs; inherit (inputs.nixpkgs) lib; };
    in
    {
      devShell.x86_64-linux = pkgs.mkShell.override { stdenv = pkgs.gcc13Stdenv; }
      {
        packages = with pkgs; [ pkg-config cmake ninja ];
        buildInputs = (with pkgs; [ eigen yaml-cpp fmt highfive tbb_2021_8.dev ])
          ++ (with localPackages; [ concurrencpp matplotplusplus zpp-bits ]);
        # hardeningDisable = [ "all" ];
        # NIX_DEBUG = "1";
      };
    };
}
