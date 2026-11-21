// @ts-check
import { defineConfig } from 'astro/config';
import sitemap from '@astrojs/sitemap';

export default defineConfig({
  site: 'https://dutra.cc',
  output: 'static',
  build: { format: 'directory' },
  // The old site had a separate publications page; send it to the section on the front page.
  redirects: { '/publications': '/#publications' },
  integrations: [sitemap({ filter: (page) => !page.includes('/publications') })],
});
