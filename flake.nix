{
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-23.05";

  outputs = inputs: let pkgs = inputs.nixpkgs.legacyPackages.x86_64-linux; in
  {
    devShell.x86_64-linux = pkgs.mkShell.override { stdenv = pkgs.gcc13Stdenv; }
    {
      buildInputs = with pkgs; [ yaml-cpp eigen fmt ];
      nativeBuildInputs = with pkgs; [ gdb ];
    };
  };
}
