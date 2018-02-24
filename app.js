img2svg = function (input, output) {
	(async () => {
		const Potrace = require("./potrace.js");

		const Jimp = require("jimp");
		const img = await Jimp.read(input);
		const bmp = img.bitmap;
		const data = Potrace.normalizePixel(bmp.data);

		//step1
		//パラメータを設定する
		Potrace.setParam({
			//Potraceパラメータ
			//NOTE:Potrace.jsではオリジナルのPotraceの全ての機能を実装していない
			//2値化する際の閾値(0.625)
			threshold: 0.5,
			//細部の省略対象とするピクセルのサイズ
			turdSize: 1,
			//パス頂点の尖りの閾値:値を0とするとカーブがなくなり, 4/3を超えると全てがカーブとして出力される.
			alphamax: 0.9,

			//Potrace.js専用・SVGパス変換用パラメータ
			//小数点下桁数(-1で無制限)
			precision: 1,
			//画像が2値化されていない場合はtrueを指定する
			noFilter: true,
			//出力結果の一部に相対パスを利用する
			isRelative: true
		});

		//step2
		//トレースを実行する
		const result = Potrace.trace(data, bmp.width, bmp.height);

		//step3
		//トレース結果からパス文字列を取得する
		const path = result.toPathString();

		//step4
		//パス文字列をSVGに埋め込んで出力する
		const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${bmp.width}" height="${bmp.height}">
<path d="${path}"/>
</svg>`;

		//output
		//process.stdout.write(svg);

		var fs = require("fs");
		fs.writeFile(output, svg);
	})();
}

for (i = 1; i < 10000; i++) {
	try {
		img2svg("D:/cut/3_png/" + i + ".png", "D:/cut/4_svg/" + i + ".svg");
	}
	catch (e) {
		print(i + "はエラー")
	}
}