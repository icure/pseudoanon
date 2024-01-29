/*
Happy linting! 💖
*/
module.exports = {
	env: {
		browser: true,
		es6: true,
		node: true,
	},
	parser: '@typescript-eslint/parser',
	parserOptions: {
		ecmaVersion: 2020, // Allows for the parsing of modern ECMAScript features
		sourceType: 'module',
	},
	plugins: ['eslint-plugin-import', 'eslint-plugin-jsdoc'],
	extends: ['plugin:@typescript-eslint/recommended', 'plugin:prettier/recommended'],
	rules: {
		'@typescript-eslint/no-explicit-any': 'off',
		'@typescript-eslint/ban-types': 'off',
	},
}
