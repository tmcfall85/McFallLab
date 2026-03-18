from pathlib import Path

from tqdm import tqdm
import subprocess

base_path = Path(
    f"/mnt/c/Users/msochor/Nextcloud/seo_tempus_fourth_delivery/vendor_tempus/"
)

# Iterate over all files and directories recursively
for path in base_path.glob("*"):
    acc_id = path.parts[-1]
    print(f"Uploading directory {acc_id}")
    inner_base_path = Path(
        f"/mnt/c/Users/msochor/Nextcloud/seo_tempus_fourth_delivery/vendor_tempus/{acc_id}/"
    )

    # Iterate over all files and directories recursively
    print(f"Uploading directory {acc_id}")
    for path in tqdm(inner_base_path.rglob("*")):
        if path.is_file():
            ssh = [
                "ssh",
                "msochor@login-hpc.rcc.mcw.edu",
                "mkdir -p /group/dseo/work/tempus_sochor_upload/"
                + "/".join(list(path.parts[6:-1])),
            ]
            scp = [
                "scp",
                str(path),
                "msochor@login-hpc.rcc.mcw.edu:/group/dseo/work/tempus_sochor_upload/"
                + "/".join(list(path.parts[6:])),
            ]
            command_parts = "\\".join(list(path.parts[6:]))
            command = f'attrib +U "C:\\Users\\msochor\\Nextcloud\\{command_parts}"'
            remove_local = [
                "/mnt/c/windows/System32/WindowsPowerShell/v1.0/powershell.exe",
                "-Command",
                command,
            ]
            subprocess.run(ssh)
            with open(path) as fp:
                try:
                    data = fp.readline()
                except:
                    data = 1
            del data
            subprocess.run(scp)
            subprocess.run(remove_local)
