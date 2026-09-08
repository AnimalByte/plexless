#![allow(dead_code)]

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;

pub struct TestDir {
    path: PathBuf,
}

impl TestDir {
    pub fn new(name: &str) -> Self {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("System clock is before UNIX epoch")
            .as_nanos();

        let path = std::env::temp_dir().join(format!(
            "simplex-test-{name}-{}-{nanos}",
            std::process::id()
        ));

        fs::create_dir_all(&path).expect("Could not create temporary test directory");

        Self { path }
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn child(&self, name: &str) -> PathBuf {
        self.path.join(name)
    }

    pub fn write(&self, name: &str, content: &str) -> PathBuf {
        let path = self.child(name);

        fs::write(&path, content).expect("Could not write test fixture");

        path
    }
}

impl Drop for TestDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}

pub fn read_gzip_text(path: &Path) -> String {
    let file = fs::File::open(path).expect("Could not open gzip output");

    let mut decoder = MultiGzDecoder::new(file);
    let mut text = String::new();

    decoder
        .read_to_string(&mut text)
        .expect("Could not decompress gzip output");

    text
}
